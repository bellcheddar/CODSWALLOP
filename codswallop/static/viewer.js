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
          // Everything in this viewport is now in the reference's coordinate frame, which is
          // the only state in which the pipeline's AlphaFold transform means anything. See
          // addAlphaFold: without this flag the model was placed in the reference's frame
          // while the viewport showed a single entry in its own, so the two came out rotated
          // against each other while both being individually correct.
          viewer._referenceFrame = entries.length ? entries[0].pdb_id : true;
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
        // The pipeline transform maps the model onto the family's REFERENCE structure, so it
        // is correct only while the viewport is in that frame. A reader who opened one entry
        // and pressed this button is looking at that entry in its own crystal frame: KRAS
        // put the model in 8T72's frame over a viewport showing 9IAY, and the two came out
        // interpenetrating and rotated. Same protein, same neighbourhood, wrong frame.
        if (af && af.u && viewer._referenceFrame) {
          fit = { u: af.u, t: af.t, source: "pipeline" };
        } else {
          // No artefact for this family, which is every family a reader assembles on the
          // live site. Fit it here instead of dropping the model in its own frame.
          fit = fitToReference(viewer, struct);
        }
        viewer._afFit = fit;
        // The state node, kept so the model can be taken off again. Without a handle on it
        // the only way to get rid of the model was to reload the whole viewer, which is why
        // the button could only ever add.
        viewer._afRef = struct.ref;
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
   * Take the AlphaFold model back off, leaving the experimental structures untouched.
   *
   * Deletes the one state node the model was built under, which takes its representation
   * and its transform with it. Returns whether there was anything to remove, so the caller
   * can keep its button label honest rather than assuming the removal happened.
   */
  function removeAlphaFold(viewer) {
    if (!viewer || !viewer._afRef) return Promise.resolve(false);
    var ref = viewer._afRef;
    viewer._afRef = null;
    viewer._afFit = null;
    return Promise.resolve()
      .then(function () { return viewer.plugin.build().delete(ref).commit(); })
      .then(function () {
        // No camera reset: the reader is looking at something, and snapping the view back
        // because a second model went away is a change they did not ask for.
        return true;
      })
      .catch(function () { return false; });   // already gone with the state tree
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

      // A predicate query, not MolScript. `lib.structure.Script` does not exist in the
      // bundled Mol*: this build exposes only structure, volume, shape, loci, math, plugin
      // and extensions, so every call here threw and the catch below returned false. The
      // Contacts panel has invited the reader to "click a residue to focus it in the 3D
      // viewer" for as long as it has existed and nothing has ever happened.
      //
      // The tests take the query CONTEXT and not a location, which is the other half of why
      // a plausible-looking rewrite still fails: read `ctx.element`.
      var S = lib.structure;
      var P = S.StructureProperties;
      var q = S.Queries.generators.atoms({
        residueTest: function (ctx) {
          return P.residue.auth_seq_id(ctx.element) === seqNumber;
        },
      });
      var sel = q(new S.QueryContext(struct));
      var loci = S.StructureSelection.toLociWithSourceUnits(sel);
      if (!loci || loci.kind === "empty-loci" || !S.StructureElement.Loci.size(loci)) {
        return false;
      }

      plugin.managers.interactivity.lociSelects.selectOnly({ loci: loci });
      plugin.managers.camera.focusLoci(loci);
      return true;
    } catch (e) {
      return false;
    }
  }

  /* Camera distance as a multiple of the scene's own radius.
     Mol*'s framing, whether through camera.reset, requestCameraReset or the toolbar's own
     "fit the visible scene into view", sits at about 2.4 radii. That is right for a protein,
     which needs room for labels and for the reader to orbit without the edges biting, and
     wrong for a fifteen-atom ligand in a 240 px drawer, where it left the molecule occupying
     a tenth of the panel.
     Set from the radius rather than as a fraction of wherever the camera currently is, so
     that applying it twice does nothing the second time. As a fraction it compounded: this
     runs on load, on two animation frames, after the drawer's transition and on every
     resize, so 0.6 became 0.6 cubed and a large ligand was clipped on all four sides while a
     small one looked right. */
  // Exactly inscribed, with no extra padding. Measured on NDG in a 433x318 drawer: at 1.08
  // the molecule filled 46% of the width, and the sphere is only about a tenth larger than
  // the molecule's own 3D diagonal, so the empty frame is the cost of bounding an elongated
  // object with a sphere rather than slack that can be tuned away. What a sphere DOES buy is
  // that nothing clips at any orientation the reader orbits to, which is the fault being
  // fixed here, so it is not traded off for a few more per cent.
  var LIGAND_MARGIN = 1.0;

  /* How far the camera has to be to hold a sphere of `radius` inside the viewport.
   *
   * Derived rather than tuned. A perspective camera sees a half-angle of fov/2, so a sphere
   * of radius r is exactly inscribed at r / sin(fov/2): at Mol*'s 45 degrees that is 2.61 r,
   * which is where its own framing sits and why that framing is right. The 1.45 r this used
   * to apply is nearer than the sphere is wide, so the molecule was pushed off every edge at
   * once. It looked like a centring fault and was a distance fault, and the two are hard to
   * tell apart when the result is a molecule running out of the frame.
   *
   * The aspect ratio is the other half. The field of view is vertical, so on a drawer panel
   * wider than it is tall the vertical direction is the tight one, and framing to the
   * horizontal half-angle clips the top and bottom of anything tall. The smaller of the two
   * is the one that has to fit.
   */
  function fitDistance(c3, radius) {
    var cam = c3.camera, st = cam.state || {};
    // Orthographic scale does not depend on distance; moving the camera would change
    // nothing except which planes clip.
    if (st.mode === "orthographic") return null;
    var fovY = st.fov || Math.PI / 4;
    var vp = cam.viewport || {};
    var aspect = (vp.width && vp.height) ? (vp.width / vp.height) : 1;
    var halfY = fovY / 2;
    var halfX = Math.atan(Math.tan(halfY) * aspect);
    var half = Math.min(halfX, halfY);
    if (!(half > 0.01)) return null;
    return (radius / Math.sin(half)) * LIGAND_MARGIN;
  }

  /** Resize the canvas to its element, then frame the scene in it. In that order. */
  function fitViewer(viewer, zoom) {
    if (!viewer) return;
    try {
      // Mol* exposes the resize differently across builds; whichever exists is the one
      // that makes the canvas agree with its container. This has to happen before the
      // camera is placed, or the projection is computed against a canvas of the wrong size
      // and the molecule sits off to one side of it.
      if (typeof viewer.handleResize === "function") viewer.handleResize();
      else if (viewer.plugin.canvas3d && viewer.plugin.canvas3d.handleResize) {
        viewer.plugin.canvas3d.handleResize();
      }
    } catch (e) { /* a build without it: the framing below is still worth doing */ }

    var c3 = viewer.plugin.canvas3d;
    // `managers.camera.reset()` restores the camera Mol* stored as its default, which is
    // computed when the scene is created and therefore before the molecule is in it: it
    // centres nothing. Reframing has to come from the scene's current bounds.
    try {
      if (c3 && typeof c3.requestCameraReset === "function") {
        c3.requestCameraReset({ moveToCenter: true });
      } else {
        viewer.plugin.managers.camera.reset();
      }
    } catch (e) { /* nothing to frame yet */ }

    if (!zoom) return;
    // A frame later, so the reset above has actually been applied: it goes through the
    // render loop rather than taking effect on the spot.
    requestAnimationFrame(function () {
      try {
        var cam = c3.camera, st = cam.state;
        var sphere = c3.boundingSphereVisible || c3.boundingSphere;
        if (!sphere || !(sphere.radius > 0)) return;
        var dx = st.position[0] - st.target[0];
        var dy = st.position[1] - st.target[1];
        var dz = st.position[2] - st.target[2];
        var now = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (!(now > 1e-6)) return;
        // Along the direction the camera is already looking from, at an absolute distance.
        var want = fitDistance(c3, sphere.radius);
        if (want == null) return;             // orthographic, or no usable viewport yet
        var k = want / now;
        cam.setState({
          position: [st.target[0] + dx * k, st.target[1] + dy * k, st.target[2] + dz * k],
        });
      } catch (e) { /* leave Mol*'s own framing */ }
    });
  }

  function c3set(viewer, props) { viewer.plugin.canvas3d.setProps(props); }

  /** Re-fit whenever the mount changes size, until it is taken off the page. */
  function watchSize(mount, onChange) {
    if (typeof ResizeObserver !== "function") return;
    var ro = new ResizeObserver(function () {
      if (!mount.isConnected) { ro.disconnect(); return; }
      onChange();
    });
    ro.observe(mount);
  }

  /**
   * One chemical component, from the RCSB's ideal-geometry file.
   *
   * Ball-and-stick rather than cartoon: a cartoon of a 30-atom ligand is nothing at all.
   * The ideal file is a few kilobytes, so this is cheap enough to open on a click, and it
   * is the real geometry rather than a 2D depiction of it.
   */
  function showLigand(host, compId) {
    host.textContent = "";
    var note = document.createElement("p");
    note.className = "vstatus";
    note.textContent = "Loading the viewer\u2026";
    host.appendChild(note);

    return ensureMolstar().then(function (molstar) {
      var mount = document.createElement("div");
      mount.className = "vmount";
      host.appendChild(mount);
      return molstar.Viewer.create(mount, {
        layoutIsExpanded: false, layoutShowControls: false, layoutShowSequence: false,
        layoutShowLog: false, layoutShowLeftPanel: false, viewportShowExpand: true,
        viewportShowSelectionMode: false, viewportShowAnimation: false,
      }).then(function (viewer) {
        var plugin = viewer.plugin;
        var url = "https://files.rcsb.org/ligands/download/" +
                  encodeURIComponent(compId.toUpperCase()) + "_ideal.sdf";
        return plugin.builders.data.download({ url: url, isBinary: false },
                                             { state: { isGhost: true } })
          .then(function (data) { return plugin.builders.structure.parseTrajectory(data, "sdf"); })
          .then(function (traj) { return plugin.builders.structure.createModel(traj); })
          .then(function (model) { return plugin.builders.structure.createStructure(model); })
          .then(function (struct) {
            return plugin.builders.structure.representation.applyPreset(struct, "atomic-detail",
                                                                        { showCarbons: true });
          })
          .then(function () {
            note.remove();
            applyTheme(viewer);
            themed.push(viewer);
            // Resize BEFORE reframing, and keep doing both while the panel is still
            // moving. `camera.reset()` frames the scene against the canvas as Mol* last
            // measured it, and the drawer this sits in slides in over 320 ms: measured at
            // the wrong moment the projection never matches the final canvas, so the
            // molecule sits off to one side and hangs over the edge. Resetting the camera
            // alone does not fix it, because the camera was never the thing that was wrong.
            // The axes gizmo is Mol*'s orientation aid for a protein. On a ligand card it
            // is a coloured artefact in the corner of a picture of one molecule.
            try {
              c3set(viewer, { camera: { helper: { axes: { name: "off", params: {} } } } });
            } catch (e) { /* a build that names it differently: leave it */ }
            fitViewer(viewer, true);
            requestAnimationFrame(function () { requestAnimationFrame(function () {
              fitViewer(viewer, true);
            }); });
            // After the drawer transition has finished, whatever its duration turned out
            // to be, and again on any later resize: expanding the viewport or rotating a
            // phone otherwise leaves the same mismatch.
            setTimeout(function () { fitViewer(viewer, true); }, 380);
            watchSize(mount, function () { fitViewer(viewer, true); });
            return viewer;
          });
      });
    });
  }

  /* ---- conservation on the surface ---------------------------------------------------
   *
   * ConSurf-style: colour the loaded structure by how conserved each residue is across the
   * whole family. The conservation is already computed per SEED position; the only hard part
   * is which residue of this structure is which position of the seed.
   *
   * That mapping is the same one that shipped wrong in the contacts artefact for 52 of 71
   * families, so it is not assumed here. Mol* gives `label_seq_id` (the entity sequence
   * index) alongside `auth_seq_id` and `label_comp_id`, and the conservation columns carry
   * the seed's own residue letter, so a candidate offset can be CHECKED: if the structure's
   * residues do not read as the seed's residues, the offset is wrong and we say so rather
   * than colouring the protein with a lie.
   */
  var THREE_TO_ONE = {
    ALA: "A", ARG: "R", ASN: "N", ASP: "D", CYS: "C", GLN: "Q", GLU: "E", GLY: "G",
    HIS: "H", ILE: "I", LEU: "L", LYS: "K", MET: "M", PHE: "F", PRO: "P", SER: "S",
    THR: "T", TRP: "W", TYR: "Y", VAL: "V", MSE: "M", SEC: "U", PYL: "O",
  };

  // Below this share of residues agreeing with the seed, the offset is not believed and
  // nothing is coloured. A real match is far above it; a wrong frame is far below.
  var MIN_AGREEMENT = 0.6;
  // How far either side of the expected offset to look. Enough for an N-terminal expression
  // tag, which is the ordinary reason a construct's numbering does not start where the
  // alignment says it does.
  var OFFSET_WINDOW = 60;

  /** Every polymer residue of the structure: chain, auth number, entity index, one letter. */
  function residueList(molstar, struct) {
    var SE = molstar.lib.structure.StructureElement;
    var P = molstar.lib.structure.StructureProperties;
    var loc = SE.Location.create(struct);
    var seen = {}, out = [];
    for (var i = 0; i < struct.units.length; i++) {
      var unit = struct.units[i];
      loc.unit = unit;
      var els = unit.elements;
      for (var j = 0; j < els.length; j++) {
        loc.element = els[j];
        if (P.atom.label_atom_id(loc) !== "CA") continue;
        var chain = P.chain.auth_asym_id(loc);
        var auth = P.residue.auth_seq_id(loc);
        var key = chain + "|" + auth;
        if (seen[key]) continue;
        seen[key] = 1;
        var one = THREE_TO_ONE[P.atom.label_comp_id(loc)];
        out.push({
          chain: chain, auth: auth,
          label: P.residue.label_seq_id(loc),
          aa: one || "X",
        });
      }
    }
    return out;
  }

  /** The offset from entity index to seed position that best explains the residues seen.
   *
   *  Scored on residue identity, not assumed from `query_beg`. An N-terminal tag shifts a
   *  construct's entity numbering away from where the alignment starts, and the expected
   *  offset is then wrong by the length of the tag. Starting the search AT the expected
   *  offset means the common case is decided on the first try and the search only earns its
   *  keep on the constructs that need it. */
  function bestOffset(residues, columns, expected) {
    var byPos = {};
    for (var i = 0; i < columns.length; i++) byPos[columns[i].pos] = columns[i].seed;
    var order = [expected], d;
    for (d = 1; d <= OFFSET_WINDOW; d++) { order.push(expected + d); order.push(expected - d); }

    var best = null;
    for (var k = 0; k < order.length; k++) {
      var off = order[k], hit = 0, tested = 0;
      for (var r = 0; r < residues.length; r++) {
        var pos = residues[r].label + off;
        var seed = byPos[pos];
        if (!seed) continue;
        tested++;
        if (seed === residues[r].aa) hit++;
      }
      if (!tested) continue;
      var score = hit / tested;
      if (!best || score > best.score) best = { offset: off, score: score, tested: tested };
      // A near-perfect frame cannot be beaten; stop rather than scanning 120 more.
      if (best.score > 0.97) break;
    }
    return best;
  }

  /** Colour the structure already loaded in `viewer` by per-seed-position conservation.
   *
   *  Resolves to a report: whether it was applied, the offset it settled on, and how well
   *  the residues agreed. The caller shows that, because a reader looking at a coloured
   *  protein deserves to know it was checked. */
  function conservationColour(viewer, pdbId, opts) {
    opts = opts || {};
    var columns = opts.columns || [];
    return ensureMolstar().then(function (molstar) {
      var plugin = viewer.plugin;
      var current = plugin.managers.structure.hierarchy.current.structures[0];
      if (!current || !current.cell || !current.cell.obj) {
        return { ok: false, reason: "no structure is loaded" };
      }
      if (!columns.length) return { ok: false, reason: "this family has no alignment" };

      var residues = residueList(molstar, current.cell.obj.data);
      if (!residues.length) return { ok: false, reason: "no polymer residues to colour" };

      var fit = bestOffset(residues, columns, (opts.queryBeg || 1) - 1);
      if (!fit || fit.score < MIN_AGREEMENT) {
        return {
          ok: false,
          reason: "this entry's residues could not be matched to the family alignment",
          agreement: fit ? fit.score : 0,
        };
      }

      var byPos = {};
      for (var i = 0; i < columns.length; i++) byPos[columns[i].pos] = columns[i];

      // Five bands rather than a continuous ramp. The underlying number is a frequency over
      // a family whose depth varies by three orders of magnitude between positions, and a
      // smooth gradient invites reading a precision into it that is not there.
      //
      // Cut at THIS FAMILY'S OWN quintiles, not at fixed thresholds. A family is a set of
      // similar sequences by construction, so conservation does not use the bottom of its
      // range: carbonic anhydrase has nothing at all below 0.5 and only 1 % below 0.7, so
      // fixed bands painted 80 % of it in the top two colours and never once used the
      // "variable" end. Quintiles guarantee the picture separates this family's own most
      // variable fifth from its most conserved, which is what the view is for. The cost is
      // that two families are not comparable by colour, so the panel says the scale is
      // relative and prints the values the bands actually fall at.
      var sorted = columns.map(function (c) { return c.conservation; })
                          .sort(function (a, b) { return a - b; });
      function quantile(q) { return sorted[Math.min(sorted.length - 1,
                                                    Math.floor(q * sorted.length))]; }
      var edges = [quantile(0.2), quantile(0.4), quantile(0.6), quantile(0.8)];
      var PALETTE = ["#4a6fa5", "#6f93bd", "#b9a06a", "#d9922f", "#e8620c"];
      var BANDS = PALETTE.map(function (colour, idx) {
        return {
          // The last band has to catch 1.0 itself, so its edge sits above the maximum.
          max: idx < edges.length ? edges[idx] : Infinity,
          colour: colour,
          label: idx === 0 ? "most variable" : (idx === 4 ? "most conserved" : ""),
        };
      });
      var buckets = BANDS.map(function () { return []; });
      var coloured = 0;
      for (var r = 0; r < residues.length; r++) {
        var col = byPos[residues[r].label + fit.offset];
        if (!col) continue;
        // Only residues that ARE the seed's residue are coloured. One that disagrees is a
        // substitution in this particular structure, and painting the family's conservation
        // onto it would be answering for a residue that is not there.
        if (col.seed !== residues[r].aa) continue;
        for (var bnd = 0; bnd < BANDS.length; bnd++) {
          if (col.conservation < BANDS[bnd].max) {
            buckets[bnd].push({ auth_asym_id: residues[r].chain, auth_seq_id: residues[r].auth });
            coloured++;
            break;
          }
        }
      }

      var colours = [{ kind: "color", params: { color: "#33405a" } }];
      BANDS.forEach(function (band, idx) {
        if (buckets[idx].length) {
          colours.push({ kind: "color",
                         params: { color: band.colour, selector: buckets[idx] } });
        }
      });

      var spec = {
        kind: "single",
        metadata: { version: "1", title: "Conservation" },
        root: { kind: "root", children: [
          { kind: "download",
            params: { url: "https://files.rcsb.org/download/" + pdbId.toUpperCase() + ".cif" },
            children: [
              { kind: "parse", params: { format: "mmcif" }, children: [
                { kind: "structure", params: { type: "model" }, children: [
                  { kind: "component", params: { selector: "polymer" }, children: [
                    { kind: "representation", params: { type: "cartoon" },
                      children: colours },
                  ] },
                  { kind: "component", params: { selector: "ligand" }, children: [
                    { kind: "representation", params: { type: "ball_and_stick" }, children: [
                      { kind: "color", params: { color: "#c9a063" } },
                    ] },
                  ] },
                ] },
              ] },
            ] },
        ] },
      };

      return viewer.loadMvsData(JSON.stringify(spec), "mvsj", { replaceExisting: true })
        .then(function () {
          return { ok: true, offset: fit.offset, agreement: fit.score,
                   coloured: coloured, residues: residues.length, bands: BANDS,
                   edges: edges };
        });
    });
  }

  global.CodswallopViewer = {
    show: show, ensure: ensureMolstar, superpose: superpose, addAlphaFold: addAlphaFold,
    removeAlphaFold: removeAlphaFold,
    focusResidue: focusResidue, showLigand: showLigand,
    conservationColour: conservationColour, _bestOffset: bestOffset,
  };
})(window);
