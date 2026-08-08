/* CODSWALLOP -- the constellation.
 *
 * A hand-rolled SVG renderer rather than a charting library, for two reasons: the node
 * pulse in the cross-highlight needs direct control of individual elements on every hover,
 * and the whole map is one radial scatter with edges, which is not worth 280 kB of D3 to
 * draw. Everything here reads the theme's CSS custom properties, so the light and dark
 * variants need no JavaScript of their own.
 *
 * Positions come from the server (layout.py), so the map, the table and any future export
 * all agree on one set of coordinates.
 */
(function (global) {
  "use strict";

  var NS = "http://www.w3.org/2000/svg";

  function el(name, attrs) {
    var node = document.createElementNS(NS, name);
    for (var k in attrs) {
      if (attrs[k] !== null && attrs[k] !== undefined) node.setAttribute(k, attrs[k]);
    }
    return node;
  }

  // X-ray cyan, cryo-EM amber, NMR violet, neutron mint. Data is cool: no brass here.
  var METHOD_VAR = {
    "X-ray": "--cyan",
    "EM": "--amber",
    "NMR": "--violet",
    "Neutron": "--mint"
  };

  function methodColour(method) {
    return "var(" + (METHOD_VAR[method] || "--mute") + ")";
  }

  /* Node radius from resolution: 1/resolution, clamped. A 0.9 A structure should read as
     bigger than a 3.5 A one without a 1.0 A outlier swallowing the field. Entries with no
     resolution (NMR, some EM) get the small end rather than vanishing.

     `scale` shrinks everything as the family grows. A 60-entry family wants generous nodes;
     drawn at that size, a 1,490-entry one is a solid carpet with no visible structure at
     all. Tied to sqrt(count) so the total ink stays roughly constant however big the family
     turns out to be. */
  function radiusFor(resolution, scale) {
    var r = resolution ? 5.2 / resolution : 2.6;
    return Math.max(1.15, Math.min(7.5, r) * (scale || 1));
  }

  function densityScale(n) {
    if (n <= 250) return 1;
    return Math.max(0.34, Math.sqrt(250 / n));
  }

  function Constellation(svg, opts) {
    this.svg = svg;
    this.opts = opts || {};
    this.nodesById = {};
    this.hot = null;
    this._size = { w: 0, h: 0 };

    var self = this;
    // The field is sized by CSS; the SVG has to follow it, and follow it again when the
    // reader turns their phone.
    this._ro = new ResizeObserver(function () { self.resize(); });
    this._ro.observe(svg.parentNode);
  }

  Constellation.prototype.render = function (map, members) {
    this.map = map;
    this.byId = {};
    members.forEach(function (m) { this.byId[m.entity_id] = m; }, this);
    this.draw();
  };

  /** How solid a node looks at depth `z`, in the -1..1 range the coordinates span. */
  function depthOpacity(z) {
    var t = Math.max(-1, Math.min(1, z));
    return (0.45 + 0.55 * (t + 1) / 2).toFixed(3);
  }

  /** A projection function for the rotation the map is at right now. */
  Constellation.prototype.projector = function () {
    var f = this._frame || { cx: 0, cy: 0, scale: 1 };
    var yaw = this.yaw || 0, pitch = this.pitch || 0;
    var cyw = Math.cos(yaw), syw = Math.sin(yaw);
    var cp = Math.cos(pitch), sp = Math.sin(pitch);
    return function (n) {
      var x = n.x, y = n.y, z = n.z || 0;
      var x1 = x * cyw + z * syw;
      var z1 = -x * syw + z * cyw;
      var y1 = y * cp - z1 * sp;
      return { x: f.cx + x1 * f.scale, y: f.cy + y1 * f.scale, z: y * sp + z1 * cp };
    };
  };

  /** Move every node to its current projection. No DOM is created or destroyed. */
  Constellation.prototype.reproject = function () {
    var self = this;
    if (!this.map || !this._frame) return;
    var project = this.projector();
    this.map.nodes.forEach(function (n) {
      var g = self.nodesById[n.id];
      if (!g) return;
      var p = project(n);
      g.setAttribute("transform", "translate(" + p.x.toFixed(2) + "," + p.y.toFixed(2) + ")");
      if (n.z != null) g.style.setProperty("--depth", depthOpacity(p.z));
    });
    // The edges are between representatives and have to follow the same rotation, or they
    // detach from the nodes they connect.
    var pos = {};
    this.map.nodes.forEach(function (n) { pos[n.id] = n; });
    var lines = this.svg.querySelectorAll(".edges line");
    for (var i = 0; i < lines.length; i++) {
      var a = pos[lines[i].getAttribute("data-a")], b = pos[lines[i].getAttribute("data-b")];
      if (!a || !b) continue;
      var pa = project(a), pb = project(b);
      lines[i].setAttribute("x1", pa.x.toFixed(2)); lines[i].setAttribute("y1", pa.y.toFixed(2));
      lines[i].setAttribute("x2", pb.x.toFixed(2)); lines[i].setAttribute("y2", pb.y.toFixed(2));
    }
  };

  Constellation.prototype.resize = function () {
    var box = this.svg.parentNode.getBoundingClientRect();
    if (!box.width || !box.height) return;
    if (Math.abs(box.width - this._size.w) < 2 && Math.abs(box.height - this._size.h) < 2) return;
    this._size = { w: box.width, h: box.height };
    if (this.map) this.draw();
  };

  Constellation.prototype.draw = function () {
    var svg = this.svg, map = this.map, self = this;
    if (!map) return;

    var box = svg.parentNode.getBoundingClientRect();
    var w = this._size.w || box.width || 600;
    var h = this._size.h || box.height || 460;
    // Reserve the strip along the bottom for the note, and keep the disc centred in what is
    // left rather than in the whole panel. Without this the largest cluster's label, printed
    // out past the rim at the six o'clock position, lands on top of the note.
    var NOTE_H = 62;
    var cx = w / 2, cy = (h - NOTE_H) / 2 + 4;

    // Fit the panel to the data rather than to the -1..1 box the coordinates are declared
    // in. A real MDS embedding almost never fills that box: a tight family occupies a
    // fraction of it and was drawn as a small knot in the middle of a large empty panel,
    // with the space going to nothing. The placeholder layout does fill it, so it is
    // unaffected, and both paths now use one rule.
    var MARGIN = 46;                       // room for a cluster label printed past the rim
    var availW = Math.max(40, w / 2 - MARGIN);
    var availH = Math.max(40, (h - NOTE_H) / 2 - MARGIN);
    // Each axis measured separately, then one isotropic scale from whichever binds. Taking
    // min(availW, availH) against the larger extent instead is over-conservative on data
    // that is not square: carbonic anhydrase spans 0.89 wide by 0.52 tall, and fitting the
    // height against the width's extent left the family using 43 % of the panel.
    var extX = 0, extY = 0;
    map.nodes.forEach(function (n) {
      extX = Math.max(extX, Math.abs(n.x || 0));
      extY = Math.max(extY, Math.abs(n.y || 0));
    });
    // A single node, or every node stacked, has no extent to fit: fall back rather than
    // dividing by zero and scaling to infinity.
    var scale = Math.min(extX > 1e-6 ? availW / extX : availW,
                         extY > 1e-6 ? availH / extY : availH);

    svg.setAttribute("viewBox", "0 0 " + w + " " + h);
    while (svg.firstChild) svg.removeChild(svg.firstChild);

    // Rotation about the vertical axis then the horizontal one, applied to the three MDS
    // coordinates before projection. Kept on the instance so a drag can re-project without
    // rebuilding the DOM: a 2,000-node family is 2,000 <g> elements, and recreating them at
    // 60 Hz is not a rotation, it is a slideshow.
    // The frame the projection needs, kept on the instance so a rotation can rebuild the
    // projector without redrawing. The trig has to be recomputed per rotation rather than
    // captured here: closing over it once meant every later projection used the angles from
    // the last full draw, so dragging updated `yaw` and moved nothing.
    self._frame = { cx: cx, cy: cy, scale: scale };
    var project = self.projector();
    var X = function (x) { return cx + x * scale; };
    var Y = function (y) { return cy + y * scale; };

    // ---- edges, behind everything -------------------------------------------------
    var pos = {};
    map.nodes.forEach(function (n) { pos[n.id] = n; });

    var gEdges = el("g", { "class": "edges" });
    map.edges.forEach(function (e) {
      var a = pos[e.a], b = pos[e.b];
      if (!a || !b) return;
      var line = el("line", {
        "class": "edge", x1: X(a.x), y1: Y(a.y), x2: X(b.x), y2: Y(b.y),
        "data-a": e.a, "data-b": e.b
      });
      gEdges.appendChild(line);
    });
    svg.appendChild(gEdges);

    // ---- cluster labels, printed on the field --------------------------------------
    // Biggest first, and a label is dropped if it would land on one already placed. A family
    // with a long tail of single-entry organisms gives those tiny clusters adjacent hairline
    // sectors, and their labels then print on top of each other; the big clusters are the
    // ones worth naming, so they get the space and the tail loses its label rather than the
    // field losing its legibility. The nodes are all still there, and hovering names them.
    var gLabels = el("g", { "class": "clusterlabels" });
    var placed = [];
    var bySize = (map.clusters || []).slice().sort(function (a, b) { return b.count - a.count; });

    bySize.forEach(function (c) {
      if (!c.count) return;
      var x = X(c.x), y = Y(c.y);
      // Keep the label inside the box: a sector pointing at the edge would otherwise print
      // its name off-canvas.
      var anchor = x < w * 0.32 ? "start" : (x > w * 0.68 ? "end" : "middle");
      var tx = Math.max(8, Math.min(w - 8, x));
      // Stop short of the note along the bottom edge, which wraps to two lines on a narrow
      // panel: a cluster label landing on top of it is the one collision this layout has.
      var ty = Math.max(14, Math.min(h - NOTE_H - 4, y));

      // Approximate the label box from its character count: measuring for real would mean
      // laying every label out and reading it back, and this only has to be good enough to
      // spot an overlap.
      var text = c.name + " · " + c.count;
      var wpx = text.length * 5.6;
      var left = anchor === "start" ? tx : (anchor === "end" ? tx - wpx : tx - wpx / 2);
      var box = { l: left, r: left + wpx, t: ty - 9, b: ty + 3 };
      var clash = placed.some(function (p) {
        return !(box.r < p.l || box.l > p.r || box.b < p.t || box.t > p.b);
      });
      if (clash) return;
      placed.push(box);

      var t = el("text", { "class": "clusterlabel", x: tx, y: ty, "text-anchor": anchor });
      t.textContent = text;
      gLabels.appendChild(t);
    });
    svg.appendChild(gLabels);

    // ---- nodes ---------------------------------------------------------------------
    var gNodes = el("g", { "class": "nodes" });
    this.nodesById = {};

    // NB: named nodeScale, not scale. `scale` above is the coordinate scale that X() and
    // Y() close over, and reassigning it here collapsed every node onto the centre point
    // while leaving the edges and labels correctly placed.
    var nodeScale = densityScale(map.nodes.length);
    // Dense families overlap heavily, so let the fill go translucent: overlap then reads as
    // density rather than as one flat block of colour.
    var fillOpacity = map.nodes.length > 600 ? 0.5 : 0.78;

    map.nodes.forEach(function (n) {
      var m = self.byId[n.id];
      if (!m) return;
      var p3 = project(n);
      var r = radiusFor(m.resolution, nodeScale);
      var g = el("g", {
        "class": "node", "data-id": n.id, tabindex: "-1",
        transform: "translate(" + p3.x.toFixed(2) + "," + p3.y.toFixed(2) + ")"
      });
      g.style.setProperty("--r0", r + "px");
      // Depth. Without it a rotation is only a reshuffle: nothing tells the eye which of
      // two overlapping nodes is in front. Opacity rather than size, because size already
      // means resolution here and must go on meaning only that.
      //
      // Through a custom property, never `style.opacity`. An inline opacity beats the
      // `.node.dim` class rule, so writing it directly silently disabled the filters: the
      // entity count changed and the map did not, which reads as the toggles being broken.
      if (n.z != null) g.style.setProperty("--depth", depthOpacity(p3.z));

      // The pulse ring: invisible until the cross-highlight fires.
      g.appendChild(el("circle", { "class": "pulse", r: r }));

      // Amber halo: this entry has a ligand somebody meant to be there. Held tight to the
      // node and thin, because in a well-studied family most entries are ligand-bound
      // (1,208 of carbonic anhydrase II's 1,490) and a generous halo turns the whole field
      // amber, which says nothing.
      if (m.has_ligand) {
        g.appendChild(el("circle", {
          "class": "halo", r: r + Math.max(1.3, 2.2 * nodeScale),
          "stroke-width": Math.max(.55, 1.1 * nodeScale)
        }));
      }

      var body = el("circle", {
        "class": "body", r: r,
        fill: methodColour(m.method), "fill-opacity": fillOpacity,
        stroke: methodColour(m.method), "stroke-width": .8 * nodeScale
      });
      g.appendChild(body);

      var title = el("title");
      title.textContent = m.pdb_id + " · " + (m.resolution ? m.resolution + " Å" : m.method)
        + " · " + (m.identity != null ? m.identity + "% id" : "")
        + (m.organism ? " · " + m.organism : "");
      g.appendChild(title);

      gNodes.appendChild(g);
      self.nodesById[n.id] = g;
    });
    svg.appendChild(gNodes);

    // ---- one delegated listener for the whole field --------------------------------
    // Rather than a listener per node: a family can have two thousand of them.
    // ---- drag to rotate --------------------------------------------------------------
    // Only when the family has a third axis. Pointer events rather than mouse, so a finger
    // works; and a drag that never moves more than a few pixels is left to fall through as
    // a click, or picking a node would become impossible.
    if (self.map.three_d) {
      var drag = null;
      svg.style.cursor = "grab";
      svg.onpointerdown = function (ev) {
        if (ev.button !== 0) return;
        // NOT preventDefault() here. Cancelling pointerdown suppresses the compatibility
        // mouse events the browser would otherwise synthesise, and `click` is one of them,
        // so it stopped every node from being openable while leaving hover working: the map
        // looked half alive. Selection is refused below on the events that actually start
        // one, which costs nothing else.
        drag = { x: ev.clientX, y: ev.clientY, moved: 0, id: ev.pointerId,
                 yaw: self.yaw || 0, pitch: self.pitch || 0 };
        // Capture keeps the drag alive when the cursor leaves the panel mid-turn. It throws
        // for a pointer id the browser has no active pointer for, so it must not be allowed
        // to take the rest of the handler down with it: without the guard the drag was
        // never registered at all.
        try { svg.setPointerCapture(ev.pointerId); } catch (e) { /* carry on uncaptured */ }
        svg.style.cursor = "grabbing";
      };
      svg.onpointermove = function (ev) {
        if (!drag) return;
        var dx = ev.clientX - drag.x, dy = ev.clientY - drag.y;
        drag.moved = Math.max(drag.moved, Math.abs(dx) + Math.abs(dy));
        self.yaw = drag.yaw + dx * 0.008;
        // Clamped, so the map cannot be turned upside down and lose its own axis labels.
        self.pitch = Math.max(-1.2, Math.min(1.2, drag.pitch + dy * 0.008));
        if (!self._raf) {
          self._raf = requestAnimationFrame(function () {
            self._raf = null;
            self.reproject();
          });
        }
      };
      var endDrag = function (ev) {
        if (!drag) return;
        var wasDrag = drag.moved > 4;
        drag = null;
        svg.style.cursor = "grab";
        try { svg.releasePointerCapture(ev.pointerId); } catch (e) { /* already released */ }
        // Swallow the click that follows a real drag, so letting go over a node does not
        // also open it.
        if (wasDrag) {
          // Swallow only the click that this drag is about to synthesise. Removed on a
          // timer as well as on use: a drag that ends outside the panel never produces one,
          // and the listener would otherwise sit there and eat the next real click.
          var once = function (e2) {
            e2.stopPropagation();
            svg.removeEventListener("click", once, true);
          };
          svg.addEventListener("click", once, true);
          setTimeout(function () { svg.removeEventListener("click", once, true); }, 350);
        }
      };
      svg.onpointerup = endDrag;
      svg.onpointercancel = endDrag;
      // The two events that begin a selection or a native image drag. Refusing these stops
      // the blue smear across the field without touching the click that follows a tap.
      svg.addEventListener("selectstart", function (ev) { ev.preventDefault(); });
      svg.addEventListener("dragstart", function (ev) { ev.preventDefault(); });
    }

    svg.onmousemove = function (ev) {
      var g = ev.target.closest ? ev.target.closest(".node") : null;
      var id = g ? g.getAttribute("data-id") : null;
      if (id !== self.hot && self.opts.onHover) self.opts.onHover(id);
    };
    svg.onmouseleave = function () { if (self.opts.onHover) self.opts.onHover(null); };
    svg.onclick = function (ev) {
      var g = ev.target.closest ? ev.target.closest(".node") : null;
      if (g && self.opts.onPick) self.opts.onPick(g.getAttribute("data-id"));
    };

    if (this._visible) this.applyVisible(this._visible);
    if (this.hot) this.setHot(this.hot);
  };

  /* Which nodes are inside the current filters. Dimming rather than removing, deliberately:
     the shape of the whole family stays legible while a filter is narrowing it. */
  Constellation.prototype.applyVisible = function (visibleSet) {
    this._visible = visibleSet;
    var ids = this.nodesById;
    for (var id in ids) {
      ids[id].classList.toggle("dim", !visibleSet.has(id));
    }
    var edges = this.svg.querySelectorAll(".edge");
    for (var i = 0; i < edges.length; i++) {
      var e = edges[i];
      var on = visibleSet.has(e.getAttribute("data-a")) && visibleSet.has(e.getAttribute("data-b"));
      e.classList.toggle("dim", !on);
    }
  };

  /* The signature: pulse one node brass. Called from the entry list, the table, and
     eventually every panel a later phase adds. */
  Constellation.prototype.setHot = function (id) {
    if (this.hot && this.nodesById[this.hot]) {
      this.nodesById[this.hot].classList.remove("hot");
    }
    this.hot = id;
    var g = id && this.nodesById[id];
    if (!g) return;
    // Re-append so the pulsing node draws above its neighbours, and so the CSS animation
    // restarts even when the same node is re-entered.
    g.parentNode.appendChild(g);
    void g.offsetWidth;
    g.classList.add("hot");
  };

  global.Constellation = Constellation;
  global.Constellation.methodColour = methodColour;
})(window);
