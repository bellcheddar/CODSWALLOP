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
import time
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
        {"id": "structures", "label": "Structures", "phase": 3, "live": True},
        {"id": "ligands", "label": "Ligands", "phase": 3, "live": True},
        {"id": "contacts", "label": "Contacts", "phase": 3, "live": True},
        {"id": "crystals", "label": "Crystals", "phase": 3, "live": True},
        {"id": "quality", "label": "Quality", "phase": 3, "live": True},
    {"id": "assembly", "label": "Assembly", "phase": 4, "live": True},
    {"id": "motifs", "label": "Motifs", "phase": 4, "live": True},
    {"id": "topology", "label": "Topology", "phase": 4, "live": True},
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
        else:
            # Stale is not the same as absent. A family that has been filed before still has
            # its row, carrying the query it was built from, so it can be rebuilt at the same
            # address without the reader having to remember what they searched for. Only a
            # slug that was never filed has nothing to go on.
            #
            # This is what makes a permalink durable. Families go stale weekly by design, and
            # a pipeline version bump stales every one of them at once, so without this a
            # shared link rots on a timer: every filed family on the live site returned "that
            # drawer is empty" the moment the parser version moved.
            known = db.family_row(slug)
            if not q and known and known.get("query"):
                q = known["query"]
            if not q:
                return render_template(
                    "message.html", heading="That drawer is empty",
                    query=slug.rsplit("-", 1)[-1],
                    message="Nothing is filed under that name yet, and a permalink alone "
                            "does not say what to rebuild it from. Search for the protein "
                            "and it will be reassembled under the same address.",
                ), 404
            display = (known or {}).get("name") or q

        return render_template("family.html", slug=slug, query=q, family=fam,
                               display_name=display)

    @app.route("/f/<slug>/dossier")
    def dossier(slug):
        """A self-contained family report: one file, no assets, prints to PDF.

        Server-rendered rather than assembled in the browser, because the point of it is to
        be a document that outlives the session it came from: something to attach to a grant
        appendix or hand to somebody starting on a target. That means no external stylesheet,
        no fonts, no scripts and no live requests, so it keeps working from a mailbox in five
        years when this app is gone.

        A family that is not cached is not rebuilt here. Assembly takes up to ninety seconds
        and this URL is the kind of thing that gets fetched by a link checker.
        """
        if not db.family_fresh(slug):
            return render_template(
                "message.html", heading="Nothing filed to report on",
                query=slug.rsplit("-", 1)[-1],
                message="A dossier is written from a family that has already been "
                        "assembled. Open the family first and the report will have "
                        "something to draw on.",
            ), 404
        fam = family_mod.decorate(db.load_family(slug))
        # Sorted here rather than in the template: a citation with no year is common enough
        # in the archive, and Jinja's sort filter compares them directly, so one undated
        # paper took the whole document down with a TypeError.
        citations = sorted(
            (c for c in (fam.get("citations") or {}).values() if c.get("title")),
            key=lambda c: (c.get("year") is None, -(c.get("year") or 0)),
        )
        return render_template("dossier.html", fam=fam, slug=slug, citations=citations,
                               generated=time.strftime("%d %B %Y", time.gmtime()))

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
        # The artefact summary rides along so a stale pipeline version is visible without a
        # shell on the droplet. `on_placeholder` is the number that matters: those families
        # render a sequence-identity placeholder where the map claims structural distance,
        # and nothing else anywhere says so.
        from . import artefacts
        try:
            art = artefacts.summary()
        except Exception:
            art = None
        return jsonify({**db.stats(), "version": config.VERSION, "artefacts": art})

    @app.route("/healthz")
    def healthz():
        return jsonify({"ok": True, "version": config.VERSION})

    return app
