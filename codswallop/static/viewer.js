/* CODSWALLOP -- the 3D viewport.
 *
 * Mol* is 5 MB, which is larger than the entire rest of the page including a 1,500-entry
 * family payload. It is therefore loaded ON DEMAND, the first time a reader actually asks
 * to see a structure, and never as part of the initial page. Somebody who came to read the
 * construct table should not pay for a renderer they never open.
 *
 * Structures are fetched by Mol* straight from the RCSB, so this needs nothing on our own
 * server: no downloads, no storage, no toolchain. That is what makes the single-entry
 * viewer separable from the superposition and interaction work, which genuinely do need a
 * local pipeline.
 */
(function (global) {
  "use strict";

  var loading = null;
  var themed = [];
  var molstarRef = null;

  /** The page's own --sky, as the 0xRRGGBB integer Mol* wants. */
  function skyColour() {
    var css = getComputedStyle(document.documentElement).getPropertyValue("--sky").trim();
    var m = /^#([0-9a-f]{6})$/i.exec(css);
    return m ? parseInt(m[1], 16) : 0x0f1522;
  }

  /* Mol* renders on a white background by default, which is a bright rectangle in the
     middle of a dark archive. There is no constructor option for it in this build, so the
     colour is set on the renderer after the plugin exists. */
  function applyTheme(viewer) {
    try {
      viewer.plugin.canvas3d.setProps({
        renderer: { backgroundColor: skyColour() },
      });
    } catch (e) { /* a Mol* build without this prop shape: leave its own default */ }
  }

  // Re-theme every open viewport when the reader flips the toggle.
  var root = document.documentElement;
  new MutationObserver(function () {
    themed.forEach(applyTheme);
  }).observe(root, { attributes: true, attributeFilter: ["data-theme"] });

  function assetUrl(name) {
    // Reuse the ?v=<mtime> the server stamped on the other assets, so the viewer is cached
    // exactly as hard as everything else and busts on the same deploy.
    var tag = document.querySelector('link[href*="theme.css"]');
    var v = tag ? (tag.getAttribute("href").split("v=")[1] || "") : "";
    return "/static/vendor/" + name + (v ? "?v=" + v : "");
  }

  /** Inject Mol* once; every later call rides the same promise. */
  function ensureMolstar() {
    if (loading) return loading;
    loading = new Promise(function (resolve, reject) {
      if (global.molstar) { molstarRef = global.molstar; return resolve(global.molstar); }
      var css = document.createElement("link");
      css.rel = "stylesheet";
      css.href = assetUrl("molstar.css");
      document.head.appendChild(css);

      var js = document.createElement("script");
      js.src = assetUrl("molstar.js");
      js.onload = function () {
        molstarRef = global.molstar;
        molstarRef ? resolve(molstarRef) : reject(new Error("Mol* did not register"));
      };
      js.onerror = function () { reject(new Error("Mol* failed to load")); };
      document.head.appendChild(js);
    });
    return loading;
  }

  /**
   * Show one PDB entry in `host`.
   * Returns a promise resolving to the viewer, so a caller can chain a colouring change.
   */
  function show(host, pdbId, opts) {
    opts = opts || {};
    host.textContent = "";
    var note = document.createElement("p");
    note.className = "vstatus";
    note.textContent = "Loading the viewer… (5 MB, once per visit)";
    host.appendChild(note);

    return ensureMolstar().then(function (molstar) {
      note.textContent = "Fetching " + pdbId + " from the RCSB…";
      var mount = document.createElement("div");
      mount.className = "vmount";
      host.appendChild(mount);

      return molstar.Viewer.create(mount, {
        layoutIsExpanded: false,
        layoutShowControls: false,
        layoutShowSequence: true,
        layoutShowLog: false,
        layoutShowLeftPanel: false,
        viewportShowExpand: true,
        viewportShowSelectionMode: false,
        viewportShowAnimation: false,
        pdbProvider: "rcsb",
        emdbProvider: "rcsb",
      }).then(function (viewer) {
        return viewer.loadPdb(pdbId).then(function () {
          note.remove();
          applyTheme(viewer);
          // The viewport must follow the page's theme toggle too, or switching to light
          // leaves a black hole in the middle of a paper-coloured page.
          themed.push(viewer);
          if (opts.onReady) opts.onReady(viewer);
          return viewer;
        });
      });
    }).catch(function (err) {
      note.className = "vstatus error";
      note.textContent = "Could not show " + pdbId + ": " + err.message +
        ". The structure is still on the RCSB: ";
      var a = document.createElement("a");
      a.href = "https://www.rcsb.org/structure/" + encodeURIComponent(pdbId);
      a.target = "_blank";
      a.rel = "noopener noreferrer";
      a.textContent = "open it there";
      note.appendChild(a);
      throw err;
    });
  }

  /* TM-align gives `u` (3x3 rotation, row-major) and `t` (translation) mapping a structure
     onto the reference. Mol* wants a COLUMN-major 4x4, so the rotation is transposed on the
     way in. Getting this backwards produces a superposition that looks plausible and is
     mirrored, which is the failure mode worth naming. */
  function mat4(u, t) {
    return [
      u[0][0], u[1][0], u[2][0], 0,
      u[0][1], u[1][1], u[2][1], 0,
      u[0][2], u[1][2], u[2][2], 0,
      t[0], t[1], t[2], 1,
    ];
  }

  /* ------------------------------------------------------------------------------------
     Superposing the AlphaFold model in the browser, when the pipeline has not done it.

     The pipeline's TM-align transform is the better answer and is used whenever it exists.
     But it only exists for a family a workstation has already embedded, and a reader who
     assembles a new family on the live site has no artefact at all: the model was being
     added in its own coordinate frame, a hundred angstroms from everything else, with a note
     admitting it. Honest, and still the wrong picture.

     This does not need TM-align. The AlphaFold model and the reference structure are the
     *same protein*, so the residue numbering is a correspondence already: match CA atoms on
     auth_seq_id and fit. No sequence alignment, no threading, nothing to get subtly wrong.
     Kabsch on those pairs is the exact least-squares rotation.
     ------------------------------------------------------------------------------------ */

  /** The Structure itself, from whatever wrapper Mol* handed back.
   *  `createStructure` resolves to a StateObjectSelector, not the data: reading `.units` off
   *  it throws, which is how this first failed. */
  function structData(x) {
    if (!x) return null;
    if (x.units) return x;
    if (x.data && x.data.units) return x.data;
    if (x.obj && x.obj.data && x.obj.data.units) return x.obj.data;
    if (x.cell && x.cell.obj && x.cell.obj.data) return x.cell.obj.data;
    return null;
  }

  /** CA coordinates by author residue number, for ONE chain of one structure.
   *
   *  One chain, not all of them. An asymmetric unit with several copies of the protein
   *  numbers each copy identically, so keying on the residue number alone silently blends
   *  coordinates from different molecules metres apart in the crystal: the fit then comes
   *  out around 2.6 A everywhere, which is too good to look broken and far worse than the
   *  sub-angstrom agreement the two structures actually have. */
  function caByResidue(molstar, struct) {
    var SE = molstar.lib.structure.StructureElement;
    var P = molstar.lib.structure.StructureProperties;
    var out = {};
    var chain = null;
    var loc = SE.Location.create(struct);
    for (var i = 0; i < struct.units.length; i++) {
      var unit = struct.units[i];
      loc.unit = unit;
      var els = unit.elements;
      for (var j = 0; j < els.length; j++) {
        loc.element = els[j];
        if (P.atom.label_atom_id(loc) !== "CA") continue;
        var asym = P.chain.auth_asym_id(loc);
        if (chain === null) chain = asym;
        if (asym !== chain) continue;
        var seq = P.residue.auth_seq_id(loc);
        if (out[seq] === undefined) {
          out[seq] = [P.atom.x(loc), P.atom.y(loc), P.atom.z(loc)];
        }
      }
    }
    return out;
  }

  /** Least-squares rotation+translation taking `mob` onto `ref` (Kabsch). */
  function kabsch(mob, ref) {
    var n = mob.length;
    if (n < 3) return null;
    var cm = [0, 0, 0], cr = [0, 0, 0], i, k;
    for (i = 0; i < n; i++) {
      for (k = 0; k < 3; k++) { cm[k] += mob[i][k] / n; cr[k] += ref[i][k] / n; }
    }
    // Covariance of the centred clouds.
    var H = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (i = 0; i < n; i++) {
      var a = [mob[i][0] - cm[0], mob[i][1] - cm[1], mob[i][2] - cm[2]];
      var b = [ref[i][0] - cr[0], ref[i][1] - cr[1], ref[i][2] - cr[2]];
      for (k = 0; k < 3; k++) {
        for (var l = 0; l < 3; l++) H[k][l] += a[k] * b[l];
      }
    }
    // Jacobi eigen-decomposition of HtH gives V; U follows from H·V·S^-1. Small and fixed
    // size, so an iterative solver is cheaper to get right than a hand-rolled SVD.
    var HtH = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (k = 0; k < 3; k++) {
      for (var m = 0; m < 3; m++) {
        var s = 0;
        for (i = 0; i < 3; i++) s += H[i][k] * H[i][m];
        HtH[k][m] = s;
      }
    }
    var ev = jacobi(HtH);
    var V = ev.vectors, w = ev.values;
    // Descending eigenvalue order, so the smallest singular value is the one we may flip.
    var idx = [0, 1, 2].sort(function (x, y) { return w[y] - w[x]; });
    var Vs = idx.map(function (c) { return [V[0][c], V[1][c], V[2][c]]; });
    var sig = idx.map(function (c) { return Math.sqrt(Math.max(w[c], 0)); });
    var Us = [];
    for (var c = 0; c < 3; c++) {
      if (sig[c] < 1e-6) { Us.push(null); continue; }
      var u = [0, 0, 0];
      for (k = 0; k < 3; k++) {
        var t = 0;
        for (m = 0; m < 3; m++) t += H[k][m] * Vs[c][m];
        u[k] = t / sig[c];
      }
      Us.push(u);
    }
    if (!Us[0] || !Us[1]) return null;
    if (!Us[2]) Us[2] = cross(Us[0], Us[1]);

    // A reflection fits the points exactly as well as a rotation and is not one. The test
    // is the determinant of the *product* U.V^T, not the handedness of U alone: checking
    // only U passed a left-handed V straight through, and this returned det -1 with an RMSD
    // of 21.8 A on a synthetic pair related by an exact 40-degree rotation.
    // R = V.U^T, which maps the mobile cloud onto the reference. Building U.V^T instead
    // gives the transpose, an equally valid-looking rotation in exactly the wrong direction.
    function build(u) {
      var M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
      for (var cc = 0; cc < 3; cc++) {
        for (var kk = 0; kk < 3; kk++) {
          for (var mm = 0; mm < 3; mm++) M[kk][mm] += Vs[cc][kk] * u[cc][mm];
        }
      }
      return M;
    }
    var R = build(Us);
    if (det3(R) < 0) {
      // Flip the column paired with the smallest singular value: it costs the least.
      Us[2] = [-Us[2][0], -Us[2][1], -Us[2][2]];
      R = build(Us);
    }
    var t2 = [0, 0, 0];
    for (k = 0; k < 3; k++) {
      t2[k] = cr[k] - (R[k][0] * cm[0] + R[k][1] * cm[1] + R[k][2] * cm[2]);
    }
    // Root-mean-square deviation after the fit, so the caller can report it rather than
    // asserting that an alignment happened.
    var sd = 0;
    for (i = 0; i < n; i++) {
      for (k = 0; k < 3; k++) {
        var p = R[k][0] * mob[i][0] + R[k][1] * mob[i][1] + R[k][2] * mob[i][2] + t2[k];
        sd += (p - ref[i][k]) * (p - ref[i][k]);
      }
    }
    return { u: R, t: t2, rmsd: Math.sqrt(sd / n), n: n };
  }

  function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
  function det3(m) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
         - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
         + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  }
  function cross(a, b) {
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  }

  /** Jacobi eigenvalue iteration for a symmetric 3x3. */
  function jacobi(A) {
    var a = [A[0].slice(), A[1].slice(), A[2].slice()];
    var v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    for (var sweep = 0; sweep < 50; sweep++) {
      var off = Math.abs(a[0][1]) + Math.abs(a[0][2]) + Math.abs(a[1][2]);
      if (off < 1e-12) break;
      for (var p = 0; p < 2; p++) {
        for (var q = p + 1; q < 3; q++) {
          if (Math.abs(a[p][q]) < 1e-14) continue;
          var theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
          var t = Math.sign(theta) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
          if (!isFinite(t)) t = 1;
          var c = 1 / Math.sqrt(t * t + 1), s = t * c;
          for (var k = 0; k < 3; k++) {
            var akp = a[k][p], akq = a[k][q];
            a[k][p] = c * akp - s * akq;
            a[k][q] = s * akp + c * akq;
          }
          for (k = 0; k < 3; k++) {
            var apk = a[p][k], aqk = a[q][k];
            a[p][k] = c * apk - s * aqk;
            a[q][k] = s * apk + c * aqk;
            var vkp = v[k][p], vkq = v[k][q];
            v[k][p] = c * vkp - s * vkq;
            v[k][q] = s * vkp + c * vkq;
          }
        }
      }
    }
    return { values: [a[0][0], a[1][1], a[2][2]], vectors: v };
  }

  /**
   * Superpose several structures in one viewport.
   * `entries` is [{pdb_id, transform:{u,t}, colour}], the first being the reference.
   */
  function superpose(host, entries, opts) {
    opts = opts || {};
    host.textContent = "";
    var note = document.createElement("p");
    note.className = "vstatus";
    note.textContent = "Loading the viewer…";
    host.appendChild(note);

    return ensureMolstar().then(function (molstar) {
      var mount = document.createElement("div");
      mount.className = "vmount";
      host.appendChild(mount);

      return molstar.Viewer.create(mount, {
        layoutIsExpanded: false, layoutShowControls: false, layoutShowSequence: false,
        layoutShowLog: false, layoutShowLeftPanel: false, viewportShowExpand: true,
        viewportShowSelectionMode: false, viewportShowAnimation: false,
        pdbProvider: "rcsb", emdbProvider: "rcsb",
      }).then(function (viewer) {
        var plugin = viewer.plugin;
        var ST = molstar.lib.plugin.StateTransforms;

        // Sequentially, not in parallel: Mol*'s state tree is a single builder and racing
        // several structure loads through it drops representations at random.
        var chain = Promise.resolve();
        entries.forEach(function (e, i) {
          chain = chain.then(function () {
            note.textContent = "Superposing " + e.pdb_id + " (" + (i + 1) + " of " +
              entries.length + ")…";
            return plugin.builders.data.download(
              { url: "https://files.rcsb.org/download/" + e.pdb_id + ".cif", isBinary: false },
              { state: { isGhost: true } }
            ).then(function (data) {
              return plugin.builders.structure.parseTrajectory(data, "mmcif");
            }).then(function (traj) {
              return plugin.builders.structure.createModel(traj);
            }).then(function (model) {
              return plugin.builders.structure.createStructure(model);
            }).then(function (struct) {
              var apply = e.transform
                ? plugin.build().to(struct).apply(ST.Model.TransformStructureConformation, {
                    transform: {
                      name: "matrix",
                      params: { data: mat4(e.transform.u, e.transform.t), transpose: false },
                    },
                  }).commit()
                : Promise.resolve();
              return apply.then(function () {
                return plugin.builders.structure.representation.applyPreset(
                  struct, "polymer-cartoon",
                  { theme: { globalName: "uniform",
                             globalColorParams: { value: e.colour } } });
              });
            });
          });
        });

        return chain.then(function () {
          note.remove();
          applyTheme(viewer);
          themed.push(viewer);
          try { plugin.managers.camera.reset(); } catch (err) { /* nothing loaded */ }
          return viewer;
        });
      });
    }).catch(function (err) {
      note.className = "vstatus error";
      note.textContent = "Could not superpose: " + err.message;
      throw err;
    });
  }

  /**
   * Add the AlphaFold model for `accession` to an existing viewer, coloured by pLDDT.
   *
   * Superposed when the pipeline supplied a transform. The alignment is done where TM-align
   * lives, in embed.py, against the same reference every experimental structure was put on,
   * so the model lands in the same frame as everything else. Without a transform (an older
   * artefact, or an accession the AFDB has no model for) it is still shown, unaligned, and
   * the caller says so rather than presenting two frames as an overlay.
   */
  /** Fit a just-loaded structure onto the first one already in the viewport.
   *  Returns the fit, or {failed: <reason>} so the caller can say why rather than shrug. */
  function fitToReference(viewer, struct) {
    try {
      var molstar = molstarRef;
      var hier = viewer.plugin.managers.structure.hierarchy.current.structures;
      if (!hier || !hier.length) return { failed: "nothing else is loaded" };
      // The reference is the first structure loaded, which `superpose` guarantees is the
      // one every pipeline transform maps onto. Anything else would fit to an arbitrary
      // member and quietly move the whole picture.
      var refData = structData(hier[0]);
      var mobData = structData(struct);
      if (!refData) return { failed: "the reference structure has no data" };
      if (!mobData) return { failed: "the model has no coordinates" };
      if (refData === mobData) return { failed: "the model is the only structure loaded" };

      var refCa = caByResidue(molstar, refData);
      var mobCa = caByResidue(molstar, mobData);

      // Residue numbering is a correspondence, but not always the same one. AlphaFold
      // numbers the full canonical sequence; a deposited entry frequently numbers the mature
      // protein, so the two are offset by the signal peptide, and pairing residue i with
      // i+1 fits everything about 2.7 A out. That is a distinctive amount of wrong: too good
      // to look broken, too poor to be the real agreement. The offset is found rather than
      // assumed, over a window wide enough for a propeptide and narrow enough that it cannot
      // wander onto a spurious register.
      var best = null;
      for (var off = -25; off <= 25; off++) {
        var m2 = [], r2 = [];
        for (var key in mobCa) {
          var partner = refCa[String(Number(key) + off)];
          if (partner) { m2.push(mobCa[key]); r2.push(partner); }
        }
        if (m2.length < 3) continue;
        var f2 = kabsch(m2, r2);
        // Most matches first, then best fit: an offset that pairs six residues perfectly is
        // not better than one that pairs two hundred well.
        if (!f2) continue;
        if (!best || m2.length > best.n * 1.15 ||
            (m2.length > best.n * 0.85 && f2.rmsd < best.rmsd)) {
          best = { fit: f2, mob: m2, ref: r2, n: m2.length, rmsd: f2.rmsd, offset: off };
        }
      }
      var mob = best ? best.mob : [], ref = best ? best.ref : [];
      var offset = best ? best.offset : 0;
      if (mob.length < 3) {
        // Author numbering that does not follow UniProt: common for a construct numbered
        // from 1, or a domain renumbered by the depositor. Named rather than shrugged at,
        // because it is a fact about the entry and not a failure of the viewer.
        return { failed: "only " + mob.length + " residue numbers are shared with " +
                         Object.keys(refCa).length + " in the structure, so the entry is " +
                         "not numbered on the canonical sequence" };
      }
      var fit = best.fit;
      if (!fit) return { failed: "the matched atoms are collinear" };

      // One round of outlier rejection. A straight fit over every matched residue is pulled
      // about by the flexible termini and by whichever loops the model placed differently:
      // carbonic anhydrase came out at 3.46 A that way, against a true core agreement well
      // under 1 A. Trimming the pairs beyond twice the RMSD and refitting reports the core,
      // which is both the honest number and the better picture.
      var keepM = [], keepR = [], cut = Math.max(2 * fit.rmsd, 1.0), i;
      for (i = 0; i < mob.length; i++) {
        var p0 = fit.u, t0 = fit.t, m = mob[i], r = ref[i], d = 0;
        for (var k = 0; k < 3; k++) {
          var v = p0[k][0] * m[0] + p0[k][1] * m[1] + p0[k][2] * m[2] + t0[k] - r[k];
          d += v * v;
        }
        if (Math.sqrt(d) <= cut) { keepM.push(m); keepR.push(r); }
      }
      if (keepM.length >= Math.max(3, 0.5 * mob.length)) {
        var refined = kabsch(keepM, keepR);
        if (refined) { refined.trimmed = mob.length - keepM.length; fit = refined; }
      }
      fit.source = "browser";
      fit.offset = offset;
      return fit;
    } catch (e) {
      // A geometry failure must not cost the reader the model itself: without a fit it
      // lands in its own frame, which is what happened before this existed.
      return { failed: "geometry error: " + e.message };
    }
  }

  function addAlphaFold(viewer, accession, af) {
    var plugin = viewer.plugin;
    // Ask the API for the file URL rather than constructing it. The model version is in the
    // filename and it moves: this was written against `-model_v4.cif` and the DB is already
    // serving v6, so every constructed URL would have 404'd.
    // The pipeline already resolved the file URL; only ask the API when it did not, which
    // is the path for a family with no embedding artefact.
    var urlPromise = (af && af.url)
      ? Promise.resolve(af.url)
      : fetch("https://alphafold.ebi.ac.uk/api/prediction/" +
              encodeURIComponent(accession.toUpperCase()))
          .then(function (r) {
            if (!r.ok) throw new Error("AlphaFold DB returned " + r.status);
            return r.json();
          })
          .then(function (rows) {
            if (!rows || !rows.length || !rows[0].cifUrl) {
              throw new Error("no model for this accession");
            }
            return rows[0].cifUrl;
          });

    var ST = molstarRef && molstarRef.lib.plugin.StateTransforms;

    return urlPromise
      .then(function (url) {
        return plugin.builders.data.download({ url: url, isBinary: false },
                                             { state: { isGhost: true } });
      })
      .then(function (data) { return plugin.builders.structure.parseTrajectory(data, "mmcif"); })
      .then(function (traj) { return plugin.builders.structure.createModel(traj); })
      .then(function (model) { return plugin.builders.structure.createStructure(model); })
      .then(function (struct) {
        // Put it on the reference before drawing it, so the camera reset below frames one
        // superposition rather than two objects a hundred angstroms apart.
        var fit = null;
        if (af && af.u) {
          fit = { u: af.u, t: af.t, source: "pipeline" };
        } else {
          // No artefact for this family, which is every family a reader assembles on the
          // live site. Fit it here instead of dropping the model in its own frame.
          fit = fitToReference(viewer, struct);
        }
        viewer._afFit = fit;
        if (fit && fit.failed) fit = null;   // nothing to apply, but the reason survives
        var placed = (fit && ST)
          ? plugin.build().to(struct).apply(ST.Model.TransformStructureConformation, {
              transform: {
                name: "matrix",
                params: { data: mat4(fit.u, fit.t), transpose: false },
              },
            }).commit()
          : Promise.resolve();
        return placed.then(function () { return struct; });
      })
      .then(function (struct) {
        // "plddt-confidence" is Mol*'s own AlphaFold theme; it reads the B-factor column,
        // which is where the AFDB stores pLDDT.
        return plugin.builders.structure.representation.applyPreset(
          struct, "polymer-cartoon",
          { theme: { globalName: "plddt-confidence" } }
        ).catch(function () {
          // A build without the confidence theme registered: fall back to a flat colour
          // rather than failing, and the caller says so in the note.
          return plugin.builders.structure.representation.applyPreset(struct, "polymer-cartoon");
        });
      })
      .then(function () {
        try { plugin.managers.camera.reset(); } catch (e) { /* nothing to frame */ }
      });
  }

  /**
   * Highlight and focus one residue (by author sequence number) in an open viewer.
   *
   * Uses Mol*'s own loci machinery rather than a representation of its own: the interaction
   * selection is transient, and adding a permanent component for every click would leave a
   * state tree full of them.
   */
  function focusResidue(viewer, seqNumber) {
    var plugin = viewer.plugin;
    var lib = molstarRef && molstarRef.lib;
    if (!lib) return false;
    try {
      var data = plugin.managers.structure.hierarchy.current.structures[0];
      if (!data || !data.cell || !data.cell.obj) return false;
      var struct = data.cell.obj.data;

      var Q = lib.structure.Script;
      var sel = Q.getStructureSelection(function (Qb) {
        return Qb.struct.generator.atomGroups({
          "residue-test": Qb.core.rel.eq([
            Qb.struct.atomProperty.macromolecular.auth_seq_id(), seqNumber,
          ]),
        });
      }, struct);
      var loci = lib.structure.StructureSelection.toLociWithSourceUnits(sel);
      if (!loci || loci.kind === "empty-loci") return false;

      plugin.managers.interactivity.lociSelects.selectOnly({ loci: loci });
      plugin.managers.camera.focusLoci(loci);
      return true;
    } catch (e) {
      return false;
    }
  }

  global.CodswallopViewer = {
    show: show, ensure: ensureMolstar, superpose: superpose, addAlphaFold: addAlphaFold,
    focusResidue: focusResidue,
  };
})(window);
