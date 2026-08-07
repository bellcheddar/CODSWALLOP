"""The Flask app.

Three request shapes, and the difference between them is the whole progressive-rendering
story:

* `/lookup?q=` resolves the input only. Cheap, cached, and the only place that can answer
  with a disambiguation card.
* `/f/<slug>` renders the drawer. If the family is already assembled it is embedded in the
  page and paints complete with no fetch at all; if not, the shell paints immediately in its
  filing state and the browser asks for the contents.
* `/api/family/<slug>` assembles and returns the family as JSON. The slow one, called only
  when there is nothing cached to serve.
"""

from __future__ import annotations

import os
from urllib.parse import urlencode

from flask import (
    Flask, abort, jsonify, redirect, render_template, request, url_for,
)

from . import config, db, family as family_mod, resolve as resolve_mod

# --------------------------------------------------------------------------------------
# The divider rail. Every phase's section from day one: the rail is the roadmap, and
# navigation never needs rebuilding as later phases ship.
# --------------------------------------------------------------------------------------
SECTIONS = [
    {
        "id": "overview", "label": "Overview", "phase": 1,
    },
    {"id": "sequences", "label": "Sequences", "phase": 2, "live": True},
    {"id": "constructs", "label": "Constructs", "phase": 2, "live": True},
    {"id": "domains", "label": "Domains", "phase": 2, "live": True},
    {
        "id": "structures", "label": "Structures", "phase": 3, "known": "methods",
        "blurb": "See it, superpose it, and make the hero map mean what it says.",
        "bullets": [
            "<b>Mol*</b> viewer, single entry and family superposition",
            "Colour by entry or by conservation; AlphaFold overlay with pLDDT colouring",
            "<b>Pairwise TM-score matrix</b> as a clustered heatmap with a dendrogram",
            "The map's node positions move to a real structural embedding, and the clusters "
            "it finds (apo vs holo, open vs closed) become first-class filters",
        ],
    },
    {
        "id": "ligands", "label": "Ligands", "phase": 3, "known": "ligands",
        "blurb": "Every chemical component in the family, and an honest verdict on which of "
                 "them anyone meant to be there.",
        "bullets": [
            "CCD ID, formula, SMILES, 2D depiction and occurrence count",
            "Explicit <b>ligand vs cryoprotectant vs buffer vs ion</b> classification: PEG, "
            "glycerol, sulphate, Tris, MES, MPD and cacodylate flagged and separable",
            "Cross-references to ChEMBL and DrugBank where they exist",
        ],
    },
    {
        "id": "contacts", "label": "Contacts", "phase": 3,
        "blurb": "PLIP run family-wide, so binding-site contacts are comparable across "
                 "entries rather than one-off.",
        "bullets": [
            "Per-entry interaction diagrams",
            "Family <b>interaction fingerprint</b> heatmap (ligand &times; residue)",
            "A &ldquo;hot residue&rdquo; ranking, mapped back onto the alignment and onto Mol*",
            "Runs as a background job with progressive fill: the map never waits for PLIP",
        ],
    },
    {
        "id": "crystals", "label": "Crystals", "phase": 3, "known": "crystals",
        "blurb": "What actually worked, parsed out of <code>_exptl_crystal_grow</code> across "
                 "the whole family.",
        "bullets": [
            "Structured precipitant, salt, buffer, pH, temperature, method and additives",
            "pH against precipitant class, coloured by resolution",
            "A &ldquo;what worked&rdquo; summary table, exported in a shape the Top96 "
            "crystallisation predictor can ingest",
        ],
    },
    {
        "id": "quality", "label": "Quality", "phase": 3, "known": "quality",
        "blurb": "A blunt traffic-light triage of which entries to trust.",
        "bullets": [
            "Clashscore, RSRZ outliers, Ramachandran and rotamer outliers",
            "R-free minus R-work gap",
            "EDS, structure-factor and raw-data availability",
        ],
    },
    {
        "id": "topology", "label": "Topology", "phase": 4,
        "blurb": "Fold cartoons generated in-house, plus everything needed to hand the "
                 "family to somebody else.",
        "bullets": [
            "Topology diagrams from DSSP as custom SVG: strands as arrows, helices as "
            "cylinders, connectivity preserved (built in-house, not scraped from PDBsum)",
            "Assembly and oligomeric state: author-assigned versus PISA-predicted, "
            "with disagreements flagged",
            "<b>Dossier export</b>: one-click self-contained HTML and PDF family report",
            "Shareable permalinks with the filter state encoded in the URL",
        ],
    },
]


def create_app() -> Flask:
    app = Flask(__name__)
    config.ensure_dirs()
    db.init()

    # ----------------------------------------------------------------------------------
    # Template helpers
    # ----------------------------------------------------------------------------------
    @app.context_processor
    def _helpers():
        def asset(filename):
            # ?v=<mtime> so a redeploy busts the browser cache even though nginx serves
            # /static/ with a one-year immutable expiry.
            try:
                v = int(os.path.getmtime(os.path.join(app.static_folder, filename)))
            except OSError:
                v = 0
            return url_for("static", filename=filename, v=v)

        return {"asset": asset, "version": config.VERSION, "sections": SECTIONS}

    @app.after_request
    def _no_stale_html(resp):
        # Flask sends no Cache-Control on a rendered template, and a browser's heuristic
        # caching then pins the old ?v= asset URLs, making a CSS or JS deploy invisible to
        # anyone who has visited before. Static files opt back in below.
        if resp.mimetype == "text/html":
            resp.headers.setdefault("Cache-Control", "no-cache, must-revalidate")
        return resp

    # ----------------------------------------------------------------------------------
    # Pages
    # ----------------------------------------------------------------------------------
    @app.route("/")
    def home():
        return render_template("landing.html", recent=db.recent_families(), query=None)

    @app.route("/lookup")
    def lookup():
        """Resolve the input, then send the reader to the family or ask which one they meant."""
        q = (request.args.get("q") or "").strip()
        if not q:
            return redirect(url_for("home"))

        try:
            result = resolve_mod.resolve(q)
        except Exception as exc:            # an upstream API is down or misbehaving
            app.logger.exception("resolve failed for %r", q)
            return render_template(
                "message.html", heading="The archive is not answering", query=q,
                message=f"Something upstream went wrong while looking that up: {exc}. "
                        f"The public APIs this runs on do occasionally wobble; try again.",
            ), 502

        if result["status"] == "ambiguous":
            return render_template("disambiguate.html", query=q, prompt=result["prompt"],
                                   candidates=result["candidates"])
        if result["status"] != "resolved":
            return render_template(
                "message.html", heading="Nothing filed under that", query=q,
                message=result.get("message", "That did not resolve to anything."),
            ), 404

        seed = result["seed"]
        slug = family_mod.slug_for(seed)
        return redirect(url_for("family_page", slug=slug, **{"q": q}))

    @app.route("/f/<slug>")
    def family_page(slug):
        """The drawer. Paints complete if the family is cached, filing if it is not."""
        q = (request.args.get("q") or "").strip()
        fam = db.load_family(slug) if db.family_fresh(slug) else None
        if fam is not None:
            fam = family_mod.decorate(fam)
            display = fam.get("name") or slug
        elif not q:
            # A cold permalink with nothing to rebuild from: the slug alone does not carry
            # enough to reconstruct the seed, so ask rather than guess at it.
            return render_template(
                "message.html", heading="That drawer is empty", query=slug.rsplit("-", 1)[-1],
                message="Nothing is filed under that name yet, and a permalink alone does not "
                        "say what to rebuild it from. Search for the protein and it will be "
                        "reassembled under the same address.",
            ), 404
        else:
            display = q

        return render_template("family.html", slug=slug, query=q, family=fam,
                               display_name=display)

    # ----------------------------------------------------------------------------------
    # API
    # ----------------------------------------------------------------------------------
    @app.route("/api/family/<slug>")
    def api_family(slug):
        """Assemble (or serve) a family as JSON. The only slow endpoint."""
        q = (request.args.get("q") or "").strip()
        force = request.args.get("refresh") == "1"

        if not force and db.family_fresh(slug):
            fam = db.load_family(slug)
            if fam is not None:
                return jsonify(family_mod.decorate(fam))

        if not q:
            abort(404, "No cached family under that name, and no query to rebuild it from.")

        try:
            result = resolve_mod.resolve(q)
            if result["status"] != "resolved":
                return jsonify({"error": result.get("message") or result.get("prompt")
                                or "That did not resolve."}), 400
            fam = family_mod.get_or_build(result["seed"], q, force=force)
        except ValueError as exc:
            return jsonify({"error": str(exc)}), 404
        except Exception as exc:
            app.logger.exception("family build failed for %r", q)
            return jsonify({"error": f"The build failed: {exc}"}), 502

        return jsonify(fam)

    @app.route("/api/stats")
    def api_stats():
        """Cache statistics. Also the per-app hit beacon for the mdeller.com launcher: it
        is a request the page's own JavaScript makes after rendering, which is what the
        launcher's nginx log regex needs to count a real page view."""
        return jsonify({**db.stats(), "version": config.VERSION})

    @app.route("/healthz")
    def healthz():
        return jsonify({"ok": True, "version": config.VERSION})

    return app
