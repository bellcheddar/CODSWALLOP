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
    var scale = Math.min(w, h - NOTE_H) / 2 - 40;

    svg.setAttribute("viewBox", "0 0 " + w + " " + h);
    while (svg.firstChild) svg.removeChild(svg.firstChild);

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
      var r = radiusFor(m.resolution, nodeScale);
      var g = el("g", {
        "class": "node", "data-id": n.id, tabindex: "-1",
        transform: "translate(" + X(n.x) + "," + Y(n.y) + ")"
      });
      g.style.setProperty("--r0", r + "px");

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
