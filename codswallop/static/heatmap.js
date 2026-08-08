/* CODSWALLOP -- the pairwise TM-score matrix, clustered.
 *
 * Clustering runs in the browser rather than in the pipeline, deliberately. The matrix is
 * already shipped in the embedding artefact, average linkage over a few hundred leaves is
 * milliseconds, and doing it here means the cut height is a live control: a reader can move
 * the threshold and watch the clusters merge, which is the only way to tell whether the
 * grouping is real or an artefact of where somebody put a line.
 *
 * Drawn on a canvas, not as SVG. A 260-representative family is 67,600 cells, and that many
 * <rect> elements makes a page a reader cannot scroll.
 */
(function (global) {
  "use strict";

  /** Average-linkage hierarchical clustering on a similarity matrix. */
  function cluster(tm) {
    var n = tm.length;
    var nodes = [];
    for (var i = 0; i < n; i++) nodes.push({ id: i, leaves: [i], height: 0, size: 1 });

    // Working distance matrix, keyed by live node index.
    var active = nodes.slice();
    var dist = [];
    for (i = 0; i < n; i++) {
      dist.push([]);
      for (var j = 0; j < n; j++) dist[i].push(1 - tm[i][j]);
    }
    var idx = {};
    active.forEach(function (a, k) { idx[a.id] = k; });

    while (active.length > 1) {
      var best = Infinity, bi = 0, bj = 1;
      for (i = 0; i < active.length; i++) {
        for (j = i + 1; j < active.length; j++) {
          var d = dist[idx[active[i].id]][idx[active[j].id]];
          if (d < best) { best = d; bi = i; bj = j; }
        }
      }
      var a = active[bi], b = active[bj];
      var merged = {
        id: n++, left: a, right: b, height: best,
        leaves: a.leaves.concat(b.leaves), size: a.size + b.size,
      };
      // Average linkage: the new row is the size-weighted mean of its children's.
      var row = [];
      for (var k = 0; k < dist.length; k++) {
        row.push((dist[idx[a.id]][k] * a.size + dist[idx[b.id]][k] * b.size) /
                 (a.size + b.size));
      }
      dist.forEach(function (r, ri) { r.push(row[ri]); });
      row.push(0);
      dist.push(row);
      idx[merged.id] = dist.length - 1;

      active.splice(bj, 1);
      active.splice(bi, 1);
      active.push(merged);
    }
    return active[0];
  }

  /** Cut the tree wherever the merge height exceeds `maxHeight`. */
  function cut(root, maxHeight) {
    var out = [];
    (function walk(node) {
      if (!node.left || node.height <= maxHeight) { out.push(node.leaves); return; }
      walk(node.left); walk(node.right);
    })(root);
    return out.sort(function (a, b) { return b.length - a.length; });
  }

  function cssVar(name) {
    return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  }

  /** Blue (dissimilar) through to mint (identical), matching the theme's data colours. */
  function colourFor(tm) {
    var t = Math.max(0, Math.min(1, (tm - 0.2) / 0.8));
    return t < 0.5
      ? mix(cssVar("--sky"), cssVar("--accent"), t * 2)
      : mix(cssVar("--accent"), cssVar("--mint"), (t - 0.5) * 2);
  }

  function mix(a, b, t) {
    function rgb(h) {
      var m = /^#?([0-9a-f]{6})$/i.exec(h.trim());
      var v = m ? parseInt(m[1], 16) : 0;
      return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
    }
    var x = rgb(a), y = rgb(b);
    return "rgb(" + x.map(function (c, i) {
      return Math.round(c + (y[i] - c) * t);
    }).join(",") + ")";
  }

  /**
   * Draw the matrix into `host`.
   * `reps` is the representative list (for labels and the cluster callback).
   */
  function draw(host, tm, reps, opts) {
    opts = opts || {};
    var n = tm.length;
    if (!n) { host.innerHTML = ""; return null; }

    var root = cluster(tm);
    var order = root.leaves;
    // Default well below the 0.5 fold boundary. Everything in a family is the same fold by
    // construction, so cutting at 0.5 returns a single cluster for almost any family; the
    // structure worth seeing (conformational states, fusion constructs) separates far
    // tighter. The beta-2 adrenergic receptor gives 1 cluster at 0.50 and 8 at 0.20.
    var height = opts.cutHeight != null ? opts.cutHeight : 0.18;
    var groups = cut(root, height);

    // Which cluster each leaf ended up in, for the callback and the colour key.
    var clusterOf = {};
    groups.forEach(function (g, gi) { g.forEach(function (leaf) { clusterOf[leaf] = gi; }); });

    var DPR = global.devicePixelRatio || 1;
    var side = Math.max(260, Math.min(620, host.clientWidth - 30));
    var cell = side / n;

    host.innerHTML = "";
    var canvas = document.createElement("canvas");
    canvas.width = side * DPR;
    canvas.height = side * DPR;
    canvas.style.width = side + "px";
    canvas.style.height = side + "px";
    canvas.className = "tmcanvas";
    host.appendChild(canvas);

    var ctx = canvas.getContext("2d");
    ctx.scale(DPR, DPR);
    for (var i = 0; i < n; i++) {
      for (var j = 0; j < n; j++) {
        ctx.fillStyle = colourFor(tm[order[i]][order[j]]);
        ctx.fillRect(j * cell, i * cell, Math.ceil(cell), Math.ceil(cell));
      }
    }
    // Cluster boundaries, so the blocks the dendrogram found are visible on the matrix.
    ctx.strokeStyle = cssVar("--brass");
    ctx.lineWidth = 1;
    var at = 0;
    groups.forEach(function (g) {
      at += g.length;
      if (at < n) {
        ctx.beginPath();
        ctx.moveTo(0, at * cell); ctx.lineTo(side, at * cell);
        ctx.moveTo(at * cell, 0); ctx.lineTo(at * cell, side);
        ctx.stroke();
      }
    });

    // Hovering a cell names the pair, which is the only way to read a 260-square matrix.
    canvas.onmousemove = function (ev) {
      var r = canvas.getBoundingClientRect();
      var ci = Math.floor((ev.clientY - r.top) / cell);
      var cj = Math.floor((ev.clientX - r.left) / cell);
      if (ci < 0 || cj < 0 || ci >= n || cj >= n) return;
      var a = reps[order[ci]], b = reps[order[cj]];
      canvas.title = a.pdb_id + " vs " + b.pdb_id + " — TM " +
        tm[order[ci]][order[cj]].toFixed(2);
      if (opts.onHover) opts.onHover(a, b);
    };
    canvas.onmouseleave = function () { if (opts.onHover) opts.onHover(null, null); };
    canvas.onclick = function (ev) {
      var r = canvas.getBoundingClientRect();
      var ci = Math.floor((ev.clientY - r.top) / cell);
      if (ci >= 0 && ci < n && opts.onPickCluster) {
        opts.onPickCluster(groups[clusterOf[order[ci]]].map(function (l) { return reps[l]; }));
      }
    };

    return { root: root, order: order, groups: groups, clusterOf: clusterOf, side: side };
  }

  global.TmHeatmap = { draw: draw, cluster: cluster, cut: cut };
})(window);
