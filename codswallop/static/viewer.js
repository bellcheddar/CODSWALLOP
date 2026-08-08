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
      if (global.molstar) return resolve(global.molstar);
      var css = document.createElement("link");
      css.rel = "stylesheet";
      css.href = assetUrl("molstar.css");
      document.head.appendChild(css);

      var js = document.createElement("script");
      js.src = assetUrl("molstar.js");
      js.onload = function () {
        global.molstar ? resolve(global.molstar) : reject(new Error("Mol* did not register"));
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
   * Not superposed onto the experimental structures: the AlphaFold DB model is in its own
   * frame, and computing the transform needs an alignment this browser has no way to do.
   * Mol* is told to reset the camera around everything, so the two sit side by side and the
   * comparison is honest about being a comparison rather than an overlay.
   */
  function addAlphaFold(viewer, accession) {
    var plugin = viewer.plugin;
    // Ask the API for the file URL rather than constructing it. The model version is in the
    // filename and it moves: this was written against `-model_v4.cif` and the DB is already
    // serving v6, so every constructed URL would have 404'd.
    return fetch("https://alphafold.ebi.ac.uk/api/prediction/" +
                 encodeURIComponent(accession.toUpperCase()))
      .then(function (r) {
        if (!r.ok) throw new Error("AlphaFold DB returned " + r.status);
        return r.json();
      })
      .then(function (rows) {
        if (!rows || !rows.length || !rows[0].cifUrl) {
          throw new Error("no model for this accession");
        }
        return plugin.builders.data.download({ url: rows[0].cifUrl, isBinary: false },
                                             { state: { isGhost: true } });
      })
      .then(function (data) { return plugin.builders.structure.parseTrajectory(data, "mmcif"); })
      .then(function (traj) { return plugin.builders.structure.createModel(traj); })
      .then(function (model) { return plugin.builders.structure.createStructure(model); })
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

  global.CodswallopViewer = {
    show: show, ensure: ensureMolstar, superpose: superpose, addAlphaFold: addAlphaFold,
  };
})(window);
