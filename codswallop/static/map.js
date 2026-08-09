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

  /* How big to draw a node that stands for `count` entries.

     Area, not radius, would be the textbook answer and is wrong here: size already means
     resolution on this map, and a stack of 486 drawn to area would be 22 times the width of
     a single structure and swallow its neighbours. The cube-rootish exponent and the hard
     cap keep a crowded point obviously crowded while leaving the resolution reading intact
     within a factor of three. */
  function stackScale(count, stackMax) {
    if (count <= 1) return 1;
    return Math.min(3.2, Math.max(1, Math.pow(count / stackMax, 0.3)));
  }

  /* Fitting the panel to the family without letting one structure decide the scale.

     A family can contain the same protein in a different fold. Hen lysozyme's 9J0L and 9J0M
     are cryo-EM amyloid fibrils: 100 % identical in sequence to the seed and TM 0.15-0.20 to
     all 76 other representatives, none of them above the 0.5 same-fold line. They sit at
     x = -1.0 while every other representative spans -0.036 to +0.046, so fitting to the
     furthest point left 1,686 of 1,688 nodes inside 5.4 % of the panel width.

     So the extent is a high quantile with generous headroom rather than the maximum, and
     anything beyond it is pinned to the rim and marked instead of setting the scale. The
     headroom is what stops this firing on an ordinary tail: it takes a point more than twice
     the 97th percentile to count as off-scale at all, and the robust extent is only adopted
     when it actually buys a materially better fit. Measured over all 68 built families:
     25 pin anything, 38 points in the whole archive, and the four lysozyme families gain
     between 12.6x and 17.2x. */
  var FIT_QUANTILE = 0.97;
  var FIT_HEADROOM = 2.0;
  var FIT_MIN_GAIN = 1.25;
  // A point has to clear the fitted extent by this much before it is called off-scale, so a
  // structure sitting essentially on the rim is not flagged as an outlier it is not.
  var FIT_OVERFLOW = 1.05;

  function robustExtent(values) {
    if (!values.length) return 0;
    var sorted = values.slice().sort(function (a, b) { return a - b; });
    var max = sorted[sorted.length - 1];
    var q = sorted[Math.floor(FIT_QUANTILE * (sorted.length - 1))] * FIT_HEADROOM;
    return (q > 0 && max / q >= FIT_MIN_GAIN) ? q : max;
  }

  function Constellation(svg, opts) {
    this.svg = svg;
    this.opts = opts || {};
    this.nodesById = {};
    // One entry per drawn element, each carrying every entity id it stands for. The filters
    // and the cross-highlight both work through this, because a drawn node is no longer
    // one-to-one with an entity once crowded points collapse.
    this.rendered = [];
    // Which stacks the reader has opened. Kept on the instance so a resize, which redraws,
    // does not silently close them again.
    this.expanded = {};
    this.hot = null;
    this._size = { w: 0, h: 0 };

    var self = this;
    // The field is sized by CSS; the SVG has to follow it, and follow it again when the
    // reader turns their phone.
    // Observes the SVG, not its container. The container also holds the field note, which
    // on a phone sits below the map and is part of that box: sizing the drawing to it made
    // the map taller by exactly the height of the caption describing it, every redraw.
    this._ro = new ResizeObserver(function () { self.resize(); });
    this._ro.observe(svg);
  }

  Constellation.prototype.render = function (map, members) {
    this.map = map;
    this.byId = {};
    // A new family is a new set of stacks: keeping the old open set would open whichever
    // points happened to share a key with the last one.
    this.expanded = {};
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

  /* Project, then hold the result inside the panel.

     Two jobs at once. A structure the fit deliberately left off the end of the scale is
     pinned to the rim rather than drawn outside the panel or allowed to set the scale; and a
     rotation, which moves points the fit never measured in that direction, can no longer
     push anything off the edge either. The vector from the centre is scaled down whole, so a
     pinned node keeps its bearing: which way it lies is the part that still means something,
     and it is the part a reader would look for. */
  Constellation.prototype.place = function (project, n) {
    var f = this._frame;
    var p = project(n);
    var dx = p.x - f.cx, dy = p.y - f.cy;
    var over = Math.max(Math.abs(dx) / f.limX, Math.abs(dy) / f.limY);
    if (over > 1) { dx /= over; dy /= over; }
    return { x: f.cx + dx, y: f.cy + dy, z: p.z };
  };

  /** Move every node to its current projection. No DOM is created or destroyed. */
  Constellation.prototype.reproject = function () {
    var self = this;
    if (!this.map || !this._frame) return;
    var project = this.projector();
    // Over the drawn elements, not over `map.nodes`: a collapsed stack has one element and
    // many nodes, so walking the nodes would move it once per member it holds and leave it
    // wherever the last of them happened to land.
    this.rendered.forEach(function (item) {
      var n = item.n;
      if (!n) return;
      var p = self.place(project, n);
      item.g.setAttribute("transform",
        "translate(" + p.x.toFixed(2) + "," + p.y.toFixed(2) + ")");
      if (n.z != null) item.g.style.setProperty("--depth", depthOpacity(p.z));
      // The off-scale chevron points outward from the centre, so a rotation that moves the
      // node has to re-aim it or it ends up pointing back into the field.
      if (item.mark) {
        var bearing = Math.atan2(p.y - self._frame.cy, p.x - self._frame.cx) * 180 / Math.PI;
        item.mark.setAttribute("transform", "rotate(" + bearing.toFixed(1)
          + ") translate(" + (item.r + 3.2).toFixed(1) + ",0)");
      }
    });
    // The cluster names ride along, or they name whatever has rotated under them.
    (this._labels || []).forEach(function (lab) {
      var p = self.place(project, lab.c);
      lab.el.setAttribute("x", p.x.toFixed(2));
      lab.el.setAttribute("y", p.y.toFixed(2));
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
    var box = this.svg.getBoundingClientRect();
    if (!box.width || !box.height) return;
    if (Math.abs(box.width - this._size.w) < 2 && Math.abs(box.height - this._size.h) < 2) return;
    this._size = { w: box.width, h: box.height };
    if (this.map) this.draw();
  };

  Constellation.prototype.draw = function () {
    var svg = this.svg, map = this.map, self = this;
    if (!map) return;

    var box = svg.getBoundingClientRect();
    var w = this._size.w || box.width || 600;
    var h = this._size.h || box.height || 460;
    // Reserve the strip along the bottom for the note, and keep the disc centred in what is
    // left rather than in the whole panel. Without this the largest cluster's label, printed
    // out past the rim at the six o'clock position, lands on top of the note.
    //
    // Measured, not assumed. This was a fixed 62px, which is two lines: any note longer than
    // that overlapped the field anyway, and on a phone the note is taken out of the overlay
    // entirely by CSS and sits below the map, where reserving anything for it wastes a
    // sixth of a short panel. Reading the element's computed position answers both without
    // the breakpoint having to be repeated here in JavaScript.
    var noteEl = svg.parentNode.querySelector(".fieldnote");
    var NOTE_H = 0;
    if (noteEl && getComputedStyle(noteEl).position === "absolute") {
      NOTE_H = Math.min(0.34 * h, noteEl.offsetHeight + 14);
    }
    var cx = w / 2, cy = (h - NOTE_H) / 2 + 4;

    /* Which points are crowded enough to draw as one node.

       Members of a family that share a construct share an MDS coordinate exactly, and the
       layout fans them around it. Past a dozen that fan is a sunflower packing whose spiral
       arms read as structure that is not there, so a crowded point is drawn once, scaled by
       how many entries it holds, and opens on click. */
    var stackMax = map.stack_max || 0;
    var byStack = {};
    if (stackMax) {
      map.nodes.forEach(function (n) {
        if (n.stack == null) return;
        (byStack[n.stack] || (byStack[n.stack] = [])).push(n);
      });
    }

    var toDraw = [];
    map.nodes.forEach(function (n) {
      var group = n.stack != null ? byStack[n.stack] : null;
      if (!group || group.length <= stackMax || self.expanded[n.stack]) {
        toDraw.push({ n: n, ids: [n.id], count: 1 });
        return;
      }
      // One member stands for the group. `si` 0 is the one the layout placed at the
      // representative's own coordinate, with no fan offset, so the stack sits where the
      // measurement actually is.
      if (n.si !== 0) return;
      toDraw.push({
        n: n, count: group.length,
        ids: group.map(function (o) { return o.id; })
      });
    });

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
    //
    // Measured over the points actually drawn rather than over every entity, so that a
    // construct solved four hundred times counts once. Weighting by entity count would let
    // one popular construct decide where the 97th percentile falls.
    var extX = robustExtent(toDraw.map(function (t) { return Math.abs(t.n.x || 0); }));
    var extY = robustExtent(toDraw.map(function (t) { return Math.abs(t.n.y || 0); }));
    // A single node, or every node stacked, has no extent to fit: fall back rather than
    // dividing by zero and scaling to infinity.
    var scale = Math.min(extX > 1e-6 ? availW / extX : availW,
                         extY > 1e-6 ? availH / extY : availH);

    // Which structures the fit has deliberately left off the end of the scale. Decided in
    // data space, not from the projected position, so that it is a property of the structure
    // and not of how far the reader has currently rotated the field.
    var offscale = 0;
    toDraw.forEach(function (t) {
      var ox = extX > 1e-6 ? Math.abs(t.n.x || 0) / extX : 0;
      var oy = extY > 1e-6 ? Math.abs(t.n.y || 0) / extY : 0;
      t.outBy = Math.max(ox, oy);
      t.offscale = t.outBy > FIT_OVERFLOW;
      if (t.offscale) offscale++;
    });
    this.offscale = offscale;

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
    self._frame = { cx: cx, cy: cy, scale: scale, limX: availW, limY: availH };
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
    // Projected, not merely scaled. On the embedded map a label names a group of nodes that
    // rotates with the rest of the field; positioning it with X()/Y() alone left the names
    // pinned to the panel while the clumps they name slid out from under them.
    self._labels = [];
    var bySize = (map.clusters || []).slice().sort(function (a, b) { return b.count - a.count; });

    bySize.forEach(function (c) {
      if (!c.count) return;
      var pl = self.place(project, c);
      var x = pl.x, y = pl.y;
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
      self._labels.push({ el: t, c: c });
    });
    svg.appendChild(gLabels);

    // ---- nodes ---------------------------------------------------------------------
    var gNodes = el("g", { "class": "nodes" });
    this.nodesById = {};


    this.rendered = [];

    // NB: named nodeScale, not scale. `scale` above is the coordinate scale that X() and
    // Y() close over, and reassigning it here collapsed every node onto the centre point
    // while leaving the edges and labels correctly placed.
    //
    // Measured against what is actually painted, not against the family's entity count.
    // Those were the same number until crowded points started collapsing; afterwards ABL1
    // was still being shrunk as though it were drawing 2,000 nodes while it drew 144, so
    // every node came out a third of the size it should be and the field looked empty.
    var nodeScale = densityScale(toDraw.length);
    // Dense families overlap heavily, so let the fill go translucent: overlap then reads as
    // density rather than as one flat block of colour.
    var fillOpacity = toDraw.length > 600 ? 0.5 : 0.78;

    // Work out how big each one is before drawing any of them, so they can be drawn largest
    // first. SVG has no z-index: the last element painted is the one on top and the one that
    // hit-testing returns. In a tight cluster four representatives can sit within seven
    // pixels of each other, and in document order a single-entry node landed squarely on the
    // centre of a 70-entry stack, so clicking the stack selected the small node instead and
    // the stack could not be opened at all. Largest at the back means whatever a reader can
    // actually see is what a click lands on.
    toDraw.forEach(function (item) {
      var m = self.byId[item.n.id];
      if (!m) { item.skip = true; return; }
      var n = item.n;
      if (item.count > 1) {
        // The stack's own reading of the family, rather than whichever entry happened to be
        // first: colour by the method most of it was solved by, halo only if most of it is
        // ligand-bound, and size from the median resolution so one 0.9 A outlier in a
        // hundred does not set the size for all of them.
        var methods = {}, withLigand = 0, resolutions = [], approxCount = 0;
        item.ids.forEach(function (id) {
          var mm = self.byId[id];
          if (!mm) return;
          methods[mm.method] = (methods[mm.method] || 0) + 1;
          if (mm.has_ligand) withLigand++;
          if (mm.resolution) resolutions.push(mm.resolution);
        });
        byStack[n.stack].forEach(function (o) { if (o.approx) approxCount++; });
        var best = null;
        for (var k in methods) { if (best === null || methods[k] > methods[best]) best = k; }
        resolutions.sort(function (a, b) { return a - b; });
        m = {
          method: best,
          resolution: resolutions.length ? resolutions[resolutions.length >> 1] : null,
          has_ligand: withLigand * 2 > item.count,
          pdb_id: m.pdb_id, organism: m.organism, identity: m.identity
        };
        n = { x: n.x, y: n.y, z: n.z, approx: approxCount * 2 > item.count };
      }
      item.m = m;
      item.node = n;
      // The radius the sort orders by has to be the radius that gets drawn, so it is derived
      // from the aggregate above rather than from the first member's resolution.
      item.r = radiusFor(m.resolution, nodeScale)
        * (item.count > 1 ? stackScale(item.count, stackMax) : 1);
    });
    toDraw = toDraw.filter(function (item) { return !item.skip; });
    toDraw.sort(function (a, b) { return b.r - a.r; });

    toDraw.forEach(function (item) {
      var n = item.node;
      var m = item.m;
      var p3 = self.place(project, n);
      var r = item.r;
      var stacked = item.count > 1;
      var g = el("g", {
        "class": "node" + (stacked ? " stacked" : "") + (n.approx ? " approx" : "")
          + (item.offscale ? " offscale" : ""),
        "data-id": item.n.id, "data-count": item.count,
        "data-stack": item.n.stack == null ? null : item.n.stack,
        tabindex: "-1",
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

      // An inferred position is drawn hollow. A third of a large family, and two thirds of
      // ABL1's, sits at a position copied from the nearest representative by identity rather
      // than measured, and the map used to draw those exactly like the ones that were
      // aligned. Set as an attribute, not from a stylesheet: a CSS rule would beat the
      // density-driven fill-opacity above and flatten every node in a crowded family.
      var body = el("circle", {
        "class": "body", r: r,
        fill: methodColour(m.method),
        "fill-opacity": n.approx ? (fillOpacity * 0.18).toFixed(3) : fillOpacity,
        stroke: methodColour(m.method), "stroke-width": .8 * nodeScale
      });
      g.appendChild(body);

      // A stack is always ringed, however few it holds, because size alone is ambiguous
      // against the resolution scale: a big node might be one good structure.
      if (stacked) {
        g.appendChild(el("circle", {
          "class": "stackring", r: r + Math.max(1.6, 2.6 * nodeScale),
          "stroke-width": Math.max(.6, 1.0 * nodeScale)
        }));
        // Past this the count is worth printing. Below it the ring is enough, and a field of
        // two-digit labels is unreadable.
        if (item.count >= 25 && r >= 5) {
          var count = el("text", { "class": "stackcount", y: r * 0.36 });
          count.textContent = item.count;
          g.appendChild(count);
        }
      }

      // A chevron pointing the way the structure actually lies, so a pinned node reads as
      // "it continues past here" rather than as a structure that genuinely sits on the rim.
      // Rotated to the bearing from the centre, which is the part of its position the fit
      // has not thrown away.
      if (item.offscale) {
        var bearing = Math.atan2(p3.y - cy, p3.x - cx) * 180 / Math.PI;
        item.mark = el("path", {
          "class": "offscalemark", d: "M0,-3.4 L3.6,0 L0,3.4",
          transform: "rotate(" + bearing.toFixed(1) + ") translate("
            + (r + 3.2).toFixed(1) + ",0)"
        });
        g.appendChild(item.mark);
      }

      var title = el("title");
      if (stacked) {
        title.textContent = item.count + " entries on one construct · "
          + (m.resolution ? "median " + m.resolution + " Å" : m.method)
          + (n.approx ? " · position inferred" : "")
          + " · click to open";
      } else {
        title.textContent = m.pdb_id + " · " + (m.resolution ? m.resolution + " Å" : m.method)
          + " · " + (m.identity != null ? m.identity + "% id" : "")
          + (m.organism ? " · " + m.organism : "")
          + (n.approx ? " · position inferred" : "");
      }
      if (item.offscale) {
        title.textContent += " · OFF SCALE: " + item.outBy.toFixed(1)
          + "x beyond the plotted range, pinned to the edge";
      }
      g.appendChild(title);

      gNodes.appendChild(g);
      // Every id the element stands for resolves to it, so hovering an entry in the table
      // still pulses the right place on the map when its point is collapsed.
      item.ids.forEach(function (id) { self.nodesById[id] = g; });
      // `n` and not `item.n`: for a stack this is the synthesised node carrying the point's
      // own coordinates, which is what a rotation has to re-project.
      self.rendered.push({ g: g, ids: item.ids, n: n, mark: item.mark, r: r });
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
        drag = { x: ev.clientX, y: ev.clientY, moved: 0, id: ev.pointerId, captured: false,
                 yaw: self.yaw || 0, pitch: self.pitch || 0 };
        // Deliberately NOT capturing here. While an element holds the pointer, the click
        // that follows is retargeted to the capturing element, so `ev.target.closest(".node")`
        // in the click handler found the <svg> and nothing was ever openable. Capture is
        // taken below, once the pointer has actually moved far enough to be a drag, which
        // is the only case that needs it.
        svg.style.cursor = "grabbing";
      };
      svg.onpointermove = function (ev) {
        if (!drag) return;
        var dx = ev.clientX - drag.x, dy = ev.clientY - drag.y;
        drag.moved = Math.max(drag.moved, Math.abs(dx) + Math.abs(dy));
        if (!drag.captured && drag.moved > 4) {
          // Now it is a drag: capture so it survives the cursor leaving the panel. It
          // throws for a pointer id with no active pointer, and must not take the handler
          // down with it.
          drag.captured = true;
          try { svg.setPointerCapture(ev.pointerId); } catch (e) { /* carry on uncaptured */ }
        }
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
        if (wasDrag) {
          try { svg.releasePointerCapture(ev.pointerId); } catch (e) { /* already released */ }
        }
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
      if (!g) {
        // Empty field closes whatever was opened. Without a way back, opening a 400-entry
        // point is a one-way trip that only a reload undoes.
        var any = false;
        for (var key in self.expanded) { any = true; break; }
        if (any) { self.expanded = {}; self.draw(); }
        return;
      }
      // A crowded point opens rather than selecting: there is no single entry to select,
      // and picking the one that happens to stand for the group would be a lie about which
      // structure the reader clicked.
      var stack = g.getAttribute("data-stack");
      if (stack && Number(g.getAttribute("data-count")) > 1) {
        self.expanded[stack] = true;
        self.draw();
        return;
      }
      if (self.opts.onPick) self.opts.onPick(g.getAttribute("data-id"));
    };

    if (this._visible) this.applyVisible(this._visible);
    if (this.hot) this.setHot(this.hot);
  };

  /* Which nodes are inside the current filters. Dimming rather than removing, deliberately:
     the shape of the whole family stays legible while a filter is narrowing it. */
  Constellation.prototype.applyVisible = function (visibleSet) {
    this._visible = visibleSet;
    // Per drawn element, and dimmed only when NONE of what it stands for survives the
    // filters. Toggling through `nodesById` instead would write the verdict of whichever
    // member came last onto the whole stack, so a point holding 400 entries would vanish
    // because one of them was filtered out.
    (this.rendered || []).forEach(function (item) {
      var on = false;
      for (var i = 0; i < item.ids.length; i++) {
        if (visibleSet.has(item.ids[i])) { on = true; break; }
      }
      item.g.classList.toggle("dim", !on);
    });
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
