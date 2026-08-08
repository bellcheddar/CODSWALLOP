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

  global.CodswallopViewer = { show: show, ensure: ensureMolstar };
})(window);
