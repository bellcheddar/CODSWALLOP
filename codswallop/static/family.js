/* CODSWALLOP -- the family page.
 *
 * ONE selection state, shared by every panel. The map, the entry list, the table, the stat
 * strip and the coverage census all read `S` and all re-render from `notify()`. Nothing
 * keeps a private copy of "what is currently showing", which is what lets a filter move in
 * one panel and land in all of them, and what every panel a later phase adds must bind to
 * rather than reinventing.
 */
(function () {
  "use strict";

  var CFG = window.CODSWALLOP;
  var $ = function (id) { return document.getElementById(id); };

  /* ==================================================================================
     1. The shared state
     ================================================================================== */
  var S = {
    family: null,
    members: [],
    visible: new Set(),      // entity ids passing the current filters
    hot: null,               // cross-highlighted entity id
    picked: null,            // entity id open in the index card
    construct: null,         // seq_id highlighted from the Constructs panel
    cluster: null,           // Set of seq_ids from a TM-matrix cluster, or null
    filters: {
      identity: 30,
      methods: null,         // Set of method names, or null for "all"
      orthologues: true,
      fusions: true,
      holoOnly: false
    },
    view: "map"
  };

  var listeners = [];
  function subscribe(fn) { listeners.push(fn); }
  function notify(what) { listeners.forEach(function (fn) { fn(what); }); }

  function recompute() {
    var f = S.filters;
    var vis = new Set();
    S.members.forEach(function (m) {
      if (m.identity != null && m.identity < f.identity) return;
      if (f.methods && !f.methods.has(methodBucket(m.method))) return;
      if (!f.orthologues && m.is_orthologue) return;
      if (!f.fusions && m.is_fusion) return;
      if (f.holoOnly && !m.has_ligand) return;
      // A cluster picked off the TM matrix filters like any other control.
      if (S.cluster && !S.cluster.has(m.seq_id)) return;
      vis.add(m.entity_id);
    });
    S.visible = vis;
  }

  /* Everything outside the four named methods is "Other": the PDB has a long tail
     (Multiple methods, fibre diffraction, solution scattering) that would otherwise put a
     dozen one-entry chips in front of the reader. */
  var METHODS = ["X-ray", "EM", "NMR", "Neutron"];
  function methodBucket(m) { return METHODS.indexOf(m) >= 0 ? m : "Other"; }

  /* ==================================================================================
     2. Small helpers
     ================================================================================== */
  function esc(s) {
    return String(s == null ? "" : s).replace(/[&<>"']/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c];
    });
  }
  function num(v, dp) {
    if (v === null || v === undefined || v === "") return "—";
    return dp === undefined ? String(v) : Number(v).toFixed(dp);
  }
  function commas(n) { return (n == null) ? "—" : n.toLocaleString("en-GB"); }
  function res(r) { return r == null ? "—" : Number(r).toFixed(2) + " Å"; }

  /* ==================================================================================
     3. Header, stat strip
     ================================================================================== */
  function renderHeader() {
    var f = S.family;
    $("famName").textContent = f.name || CFG.slug;
    $("stampCount").textContent = commas(f.stats.entries) + " ENTRIES";
    document.title = (f.name || CFG.slug) + " — CODSWALLOP";

    var bits = [];
    if (f.organism) bits.push('<span class="tag">' + esc(f.organism) + "</span>");
    if (f.seed) bits.push('<span class="tag">seed ' + esc(f.seed) + "</span>");
    (f.pfam || []).slice(0, 3).forEach(function (p) {
      bits.push('<a class="tag" href="https://www.ebi.ac.uk/interpro/entry/pfam/' +
        encodeURIComponent(p.id) + '/" target="_blank" rel="noopener noreferrer">' +
        esc(p.id) + " " + esc(p.name || "") + "</a>");
    });
    (f.interpro || []).slice(0, 2).forEach(function (p) {
      bits.push('<a class="tag" href="https://www.ebi.ac.uk/interpro/entry/InterPro/' +
        encodeURIComponent(p.id) + '/" target="_blank" rel="noopener noreferrer">' +
        esc(p.id) + "</a>");
    });
    bits.push('<span class="tag">' + commas(f.seed_length) + " aa seed</span>");
    $("famMeta").innerHTML = bits.join("");

    var note = f.note || "";
    if (f.truncated) {
      note += " The PDB holds " + commas(f.total_hits) + " matching entities; the closest " +
        commas(f.stats.entities) + " were filed.";
    }
    $("famNote").textContent = note;
  }

  function renderStats() {
    var s = S.family.stats, c = s.coverage;
    var tiles = [
      { k: "Entries", v: commas(s.entries) },
      { k: "Polymer entities", v: commas(s.entities) },
      { k: "Distinct constructs", v: commas(s.constructs) },
      { k: "Organisms", v: commas(s.organisms) },
      { k: "Ligand-bound", v: commas(s.holo_entries) },
      { k: "Distinct ligands", v: commas(s.ligands) },
      { k: "Best resolution", v: res(s.best_resolution) },
      { k: "Median resolution", v: res(s.median_resolution) },
      // Prefer the density figure once Phase 2 has measured it: "present in the construct
      // but not in the density" is the question a construct designer actually has, and the
      // construct-level figure reads 0 % for any well-studied protein. Falls back to the
      // Phase 1 figure, labelled differently, when no per-chain data was available.
      (c.rarely_resolved_pct != null ? {
        k: "Rarely resolved", flag: true,
        v: c.rarely_resolved_pct + '%<small> of seed</small>',
        title: c.rarely_resolved + " of " + c.length + " seed residues are resolved in under " +
          Math.round(c.rarely_cut * 100) + "% of the constructs that contained them. " +
          "Present in the crystal, absent from the density."
      } : {
        k: "Thinly covered", flag: true,
        v: c.thin_pct + '%<small> of seed</small>',
        title: "Seed residues present in fewer than " + c.thin_cut + " of the family's " +
          "constructs (" + c.thin + " of " + c.length + " residues). Per-chain density was " +
          "not available for this family."
      }),
      {
        k: "Median construct", v: c.median_coverage + '%<small> of seed</small>',
        title: "The median deposited construct contains this fraction of the seed sequence."
      }
    ];
    $("statStrip").innerHTML = tiles.map(function (t) {
      return '<div class="stat' + (t.flag ? " flag" : "") + '"' +
        (t.title ? ' title="' + esc(t.title) + '"' : "") +
        '><div class="v">' + t.v + '</div><div class="k">' + esc(t.k) + "</div></div>";
    }).join("");
  }

  /* ==================================================================================
     4. Filters
     ================================================================================== */
  function renderFilters() {
    var s = S.family.stats;

    // The slider bounds itself to the identity range the family actually spans, so it never
    // offers a range that filters out everything or nothing across most of its travel.
    var lo = Math.floor(s.identity_min == null ? 30 : s.identity_min);
    var hi = Math.ceil(s.identity_max == null ? 100 : s.identity_max);
    var slider = $("identity");
    slider.min = lo;
    slider.max = Math.max(hi, lo + 1);
    slider.value = lo;
    S.filters.identity = lo;
    $("identityVal").textContent = "≥" + lo + "%";

    // Method chips, from what the family actually contains.
    var counts = {};
    S.members.forEach(function (m) {
      var b = methodBucket(m.method);
      counts[b] = (counts[b] || 0) + 1;
    });
    var box = $("methodChips");
    box.querySelectorAll(".chip").forEach(function (n) { n.remove(); });
    Object.keys(counts).sort(function (a, b) { return counts[b] - counts[a]; })
      .forEach(function (name) {
        var b = document.createElement("button");
        b.type = "button";
        b.className = "chip";
        b.dataset.method = name;
        b.setAttribute("aria-pressed", "true");
        b.innerHTML = '<span class="dot"></span>' + esc(name) +
          ' <span class="n">' + commas(counts[name]) + "</span>";
        b.addEventListener("click", function () { toggleMethod(name, b); });
        box.appendChild(b);
      });

    $("nOrth").textContent = "(" + commas(s.orthologues) + ")";
    $("nFusion").textContent = "(" + commas(s.fusions) + ")";
    // Nothing to include or exclude: say so by disabling rather than offering a control
    // that cannot change the answer.
    $("optOrth").disabled = !s.orthologues;
    $("optFusion").disabled = !s.fusions;
  }

  function toggleMethod(name, btn) {
    var all = Array.prototype.map.call(
      document.querySelectorAll("#methodChips .chip"), function (c) { return c.dataset.method; });
    if (!S.filters.methods) S.filters.methods = new Set(all);
    if (S.filters.methods.has(name)) S.filters.methods.delete(name);
    else S.filters.methods.add(name);
    // Empty means nothing shows, which is never what a click on the last chip meant.
    if (S.filters.methods.size === 0) S.filters.methods = new Set(all);
    document.querySelectorAll("#methodChips .chip").forEach(function (c) {
      c.setAttribute("aria-pressed", S.filters.methods.has(c.dataset.method) ? "true" : "false");
    });
    apply();
  }

  function wireFilters() {
    $("identity").addEventListener("input", function () {
      S.filters.identity = Number(this.value);
      $("identityVal").textContent = "≥" + this.value + "%";
      apply();
    });
    $("optOrth").addEventListener("change", function () {
      S.filters.orthologues = this.checked; apply();
    });
    $("optFusion").addEventListener("change", function () {
      S.filters.fusions = this.checked; apply();
    });
    $("optHolo").addEventListener("change", function () {
      S.filters.holoOnly = this.checked; apply();
    });
  }

  function apply() {
    recompute();
    notify("visible");
  }

  /* ==================================================================================
     5. The entry list, permanently beside the map
     ================================================================================== */
  var rowsById = {};

  function renderList() {
    var box = $("entryRows");
    box.innerHTML = "";
    rowsById = {};

    var shown = S.members.filter(function (m) { return S.visible.has(m.entity_id); });
    $("listCount").textContent = commas(shown.length);

    if (!shown.length) {
      box.innerHTML = '<p class="empty">Nothing passes those filters. ' +
        "Widen the identity threshold or put a method back.</p>";
      return;
    }

    // Cap what is in the DOM. Two thousand rows is slower to scroll than it is useful, and
    // the table view exists for anyone who wants every row at once.
    var LIMIT = 300;
    var frag = document.createDocumentFragment();
    shown.slice(0, LIMIT).forEach(function (m) {
      var row = document.createElement("div");
      row.className = "erow" + (m.has_ligand ? " holo" : "");
      row.setAttribute("data-id", m.entity_id);
      row.setAttribute("tabindex", "0");
      row.setAttribute("role", "button");
      row.innerHTML =
        '<span class="pdbid"><span class="mdot" style="background:' +
          Constellation.methodColour(m.method) + '"></span>' + esc(m.pdb_id) + "</span>" +
        '<span class="desc">' + esc(m.description || m.title || "") + "</span>" +
        '<span class="res">' + (m.resolution ? Number(m.resolution).toFixed(2) : "—") + "</span>";
      frag.appendChild(row);
      rowsById[m.entity_id] = row;
    });
    box.appendChild(frag);

    if (shown.length > LIMIT) {
      var more = document.createElement("p");
      more.className = "empty";
      more.textContent = commas(shown.length - LIMIT) + " more — open the table view for all of them.";
      box.appendChild(more);
    }
  }

  function wireList() {
    var box = $("entryRows");
    // Delegated: the list is re-rendered on every filter change.
    box.addEventListener("mouseover", function (ev) {
      var row = ev.target.closest(".erow");
      if (row) setHot(row.getAttribute("data-id"));
    });
    box.addEventListener("mouseleave", function () { setHot(null); });
    box.addEventListener("focusin", function (ev) {
      var row = ev.target.closest(".erow");
      if (row) setHot(row.getAttribute("data-id"));
    });
    box.addEventListener("click", function (ev) {
      var row = ev.target.closest(".erow");
      if (row) openCard(row.getAttribute("data-id"));
    });
    box.addEventListener("keydown", function (ev) {
      var row = ev.target.closest(".erow");
      if (row && (ev.key === "Enter" || ev.key === " ")) {
        ev.preventDefault();
        openCard(row.getAttribute("data-id"));
      }
    });
  }

  /* ==================================================================================
     6. The cross-highlight. One state, two representations.
     ================================================================================== */
  var constellation = null;

  function setHot(id) {
    if (S.hot === id) return;
    if (S.hot && rowsById[S.hot]) rowsById[S.hot].classList.remove("hot");
    S.hot = id;
    if (id && rowsById[id]) {
      rowsById[id].classList.add("hot");
      // Only scroll when the row is actually out of view: a list that jumps under the
      // cursor on every hover is worse than one that occasionally does not follow.
      var r = rowsById[id].getBoundingClientRect();
      var b = $("entryRows").getBoundingClientRect();
      if (r.top < b.top || r.bottom > b.bottom) {
        rowsById[id].scrollIntoView({ block: "nearest" });
      }
    }
    if (constellation) constellation.setHot(id);
    highlightTableRow(id);
  }

  /* ==================================================================================
     7. The index card
     ================================================================================== */
  function memberById(id) {
    for (var i = 0; i < S.members.length; i++) {
      if (S.members[i].entity_id === id) return S.members[i];
    }
    return null;
  }

  function openCard(id) {
    var m = memberById(id);
    if (!m) return;
    S.picked = id;

    $("cardTitle").textContent = m.pdb_id + " · entity " + m.entity_id.split("_")[1];
    $("cardTitleText").textContent = m.title || "";

    var cell = m.cell;
    var rows = [
      ["Description", m.description, true],
      ["Method", m.method],
      ["Resolution", m.resolution ? Number(m.resolution).toFixed(2) + " Å" : "—"],
      ["R-work / R-free", (m.r_work != null ? Number(m.r_work).toFixed(3) : "—") + " / " +
        (m.r_free != null ? Number(m.r_free).toFixed(3) : "—")],
      ["Space group", m.space_group || "—"],
      ["Unit cell", cell ? [cell.a, cell.b, cell.c].map(function (v) {
        return v == null ? "?" : Number(v).toFixed(1);
      }).join(" × ") + " Å, " + [cell.alpha, cell.beta, cell.gamma].map(function (v) {
        return v == null ? "?" : Number(v).toFixed(1);
      }).join(" / ") + "°" : "—"],
      ["Identity to seed", m.identity != null ? m.identity + "%" : "—"],
      ["Length", m.seq_length ? m.seq_length + " aa" : "—"],
      ["Chains", (m.chains || []).join(", ") || "—"],
      ["Assemblies", m.assembly_count],
      ["Source", m.organism, true],
      ["Expressed in", m.host_organism || "—", true],
      ["UniProt", m.uniprot || "—"],
      ["Deposited", m.deposit_date || "—"],
      ["Released", m.release_date || "—"]
    ];

    var html = '<dl class="kv">' + rows.map(function (r) {
      return "<dt>" + esc(r[0]) + '</dt><dd class="' + (r[2] ? "prose" : "") + '">' +
        esc(r[1] == null || r[1] === "" ? "—" : r[1]) + "</dd>";
    }).join("") + "</dl>";

    // Construct flags: the Phase 1 version of what the Constructs panel will say properly.
    var flags = [];
    if (m.is_fusion) {
      flags.push("Carries an apparent fusion partner or extension: the deposited construct " +
        "is " + m.seq_length + " residues, of which " + m.aligned_length + " align to the seed.");
    }
    if (m.is_orthologue) flags.push("Orthologue: a different source organism from the seed.");
    if (flags.length) {
      html += '<div class="csect"><h3>Construct notes</h3><p style="color:var(--dim);' +
        'font-size:12.5px;line-height:1.6;margin:0">' +
        flags.map(esc).join("<br>") + "</p></div>";
    }

    if ((m.ligands || []).length) {
      html += '<div class="csect"><h3>Chemical components (' + m.ligands.length +
        ')</h3><div class="liglist">' + m.ligands.map(function (id) {
          var comp = S.family.components[id] || {};
          var real = NON_LIGANDS.indexOf(id) < 0;
          return '<span class="lig' + (real ? " real" : "") + '" title="' +
            esc(comp.name || "") + '">' + esc(id) + "</span>";
        }).join("") + "</div>" +
        '<p style="color:var(--mute);font-size:11px;margin:8px 0 0">' +
        "Amber components are not on the buffer, ion and cryoprotectant list. " +
        "Phase 3 classifies every one of them properly.</p></div>";
    }

    var c = m.cite_id ? S.family.citations[m.cite_id] : null;
    if (c) {
      html += '<div class="csect"><h3>Primary citation</h3>' +
        '<p style="color:var(--dim);font-size:12.5px;line-height:1.6;margin:0">' +
        esc(c.title || "") + "<br><i>" + esc(c.journal || "") + "</i> " +
        esc(c.volume || "") + " " + esc(c.pages || "") + " (" + esc(c.year || "") + ")" +
        (c.authors && c.authors.length ? "<br>" + esc(c.authors.slice(0, 6).join(", ")) +
          (c.authors.length > 6 ? " et al." : "") : "") +
        (c.doi ? '<br><a href="https://doi.org/' + encodeURIComponent(c.doi) +
          '" target="_blank" rel="noopener noreferrer">' + esc(c.doi) + "</a>" : "") +
        "</p></div>";
    }

    var seq = m.seq_id ? S.family.sequences[m.seq_id] : null;
    if (seq) {
      html += '<div class="csect"><h3>Deposited sequence (' + seq.length +
        ' aa)</h3><div class="seqbox">' + esc(seq) + "</div></div>";
    }

    html += '<div class="csect"><h3>Structure</h3>' +
      '<button class="btn ghost" type="button" id="cardShow3D">Show ' + esc(m.pdb_id) +
      ' in 3D</button>' +
      '<div class="vhost card" id="cardViewer" hidden></div></div>';

    html += '<div class="outlinks">' +
      '<a href="https://www.rcsb.org/structure/' + encodeURIComponent(m.pdb_id) +
        '" target="_blank" rel="noopener noreferrer">RCSB</a>' +
      '<a href="https://www.ebi.ac.uk/pdbe/entry/pdb/' + encodeURIComponent(m.pdb_id.toLowerCase()) +
        '" target="_blank" rel="noopener noreferrer">PDBe</a>' +
      (m.uniprot ? '<a href="https://www.uniprot.org/uniprotkb/' + encodeURIComponent(m.uniprot) +
        '" target="_blank" rel="noopener noreferrer">UniProt</a>' : "") +
      '<a href="https://files.rcsb.org/download/' + encodeURIComponent(m.pdb_id) +
        '.cif" target="_blank" rel="noopener noreferrer">mmCIF</a>' +
      "</div>";

    $("cardBody").innerHTML = html;
    // On demand, not on open: a reader flicking through entries should not pull 5 MB per card.
    var show3d = $("cardShow3D");
    if (show3d) {
      show3d.addEventListener("click", function () {
        var host = $("cardViewer");
        host.hidden = false;
        show3d.disabled = true;
        window.CodswallopViewer.show(host, m.pdb_id).catch(function () {
          show3d.disabled = false;
        });
      });
    }
    $("indexCard").hidden = false;
    // One frame before adding the class, so the transform transition actually runs rather
    // than being skipped as part of the same style recalculation that unhides the panel.
    requestAnimationFrame(function () {
      $("indexCard").classList.add("on");
      $("scrim").classList.add("on");
    });
    $("cardClose").focus();
  }

  function closeCard() {
    $("indexCard").classList.remove("on");
    $("scrim").classList.remove("on");
    S.picked = null;
    setTimeout(function () { $("indexCard").hidden = true; }, 320);
  }

  // Kept in step with _NON_LIGANDS in family.py; used only to tint the card's chips.
  var NON_LIGANDS = ["HOH", "DOD", "SO4", "PO4", "CL", "NA", "K", "MG", "CA", "ZN", "MN",
    "FE", "NI", "CO", "CU", "CD", "HG", "BR", "IOD", "ACT", "EDO", "GOL", "PEG", "PG4",
    "PGE", "1PE", "2PE", "MPD", "DMS", "TRS", "MES", "EPE", "IMD", "FMT", "CIT", "TAR",
    "MLA", "ACY", "NH4", "CAC", "BME", "DTT", "TCE", "AZI", "NO3", "CO3", "F", "LI", "CS",
    "RB", "SR", "BA"];

  /* ==================================================================================
     8. Coverage census
     ================================================================================== */
  function renderCoverage() {
    var c = S.family.stats.coverage;
    var svg = $("coverage");
    var depth = c.depth || [];
    if (!depth.length) { svg.innerHTML = ""; return; }

    var W = 1000, H = 54, pad = 14;
    var n = depth.length, maxd = c.max_depth || 1;
    svg.setAttribute("viewBox", "0 0 " + W + " " + H);
    svg.setAttribute("preserveAspectRatio", "none");

    var bw = (W - pad * 2) / n;
    var parts = [];
    for (var i = 0; i < n; i++) {
      var d = depth[i];
      var h = Math.max(d > 0 ? 1 : 0, ((H - 16) * d) / maxd);
      var cls = d === 0 ? "none" : (d < c.thin_cut ? "thin" : "depth");
      // A zero-depth residue still needs to be visible, so it draws a full-height ghost.
      if (d === 0) h = H - 16;
      parts.push('<rect class="' + cls + '" x="' + (pad + i * bw).toFixed(2) +
        '" y="' + (H - 14 - h).toFixed(2) + '" width="' + Math.max(bw, .6).toFixed(2) +
        '" height="' + h.toFixed(2) + '"><title>residue ' + (i + 1) + ": " + d +
        " constructs</title></rect>");
    }
    parts.push('<text class="covaxis" x="' + pad + '" y="' + (H - 3) + '">1</text>');
    parts.push('<text class="covaxis" x="' + (W - pad) + '" y="' + (H - 3) +
      '" text-anchor="end">' + n + "</text>");
    svg.innerHTML = parts.join("");

    $("covLegend").innerHTML =
      '<span><i style="background:var(--accent)"></i>in ≥' + c.thin_cut + " constructs</span>" +
      '<span><i style="background:var(--oxide)"></i>thin: fewer than ' + c.thin_cut +
        " (" + c.thin + " residues)</span>" +
      '<span><i style="background:var(--oxide);opacity:.45"></i>in none (' + c.uncovered + ")</span>" +
      "<span>Median construct covers " + c.median_coverage + "% of the seed.</span>" +
      (c.gaps && c.gaps.length ? "<span>Poorest region: " + c.gaps[0].start + "–" +
        c.gaps[0].end + " (" + c.gaps[0].depth + " constructs).</span>" : "");
  }

  /* ==================================================================================
     9. The table (Tabulator), sortable, filterable, exportable
     ================================================================================== */
  var table = null;

  function tableData() {
    return S.members.filter(function (m) { return S.visible.has(m.entity_id); })
      .map(function (m) {
        return {
          entity_id: m.entity_id, pdb_id: m.pdb_id,
          description: m.description || m.title || "",
          method: m.method, resolution: m.resolution,
          identity: m.identity, r_work: m.r_work, r_free: m.r_free,
          space_group: m.space_group, length: m.seq_length,
          chains: (m.chains || []).join(" "), organism: m.organism,
          host: m.host_organism, uniprot: m.uniprot,
          ligands: (m.ligands || []).join(" "),
          holo: m.has_ligand ? "yes" : "no",
          fusion: m.is_fusion ? "yes" : "no",
          deposited: m.deposit_date, released: m.release_date,
          title: m.title || ""
        };
      });
  }

  function buildTable() {
    if (table) { table.replaceData(tableData()); return; }

    var mono = function (cell) { cell.getElement().classList.add("mono"); return cell.getValue(); };
    // Fixed decimals, because the archive reports them at whatever precision the depositor
    // felt like: an R-work of 0.11157 beside one of 0.16 is noise, not extra information.
    var fixed = function (dp) {
      return function (cell) {
        cell.getElement().classList.add("mono");
        var v = cell.getValue();
        return (v === null || v === undefined || v === "") ? "" : Number(v).toFixed(dp);
      };
    };

    table = new Tabulator("#entryTable", {
      data: tableData(),
      layout: "fitDataFill",
      height: "620px",
      index: "entity_id",
      placeholder: "Nothing passes those filters.",
      columnDefaults: { headerFilterLiveFilter: true, resizable: true },
      columns: [
        {
          title: "PDB", field: "pdb_id", width: 78, headerFilter: "input", frozen: true,
          formatter: function (cell) {
            cell.getElement().classList.add("idcell");
            return cell.getValue();
          }
        },
        { title: "Entity", field: "entity_id", width: 92, visible: false, formatter: mono },
        { title: "Description", field: "description", width: 260, headerFilter: "input" },
        { title: "Method", field: "method", width: 104, headerFilter: "list",
          headerFilterParams: { valuesLookup: true, clearable: true } },
        { title: "Res (Å)", field: "resolution", width: 88, hozAlign: "right", formatter: fixed(2),
          sorter: "number" },
        { title: "Id (%)", field: "identity", width: 82, hozAlign: "right", formatter: fixed(1),
          sorter: "number" },
        { title: "R-work", field: "r_work", width: 84, hozAlign: "right", formatter: fixed(3),
          sorter: "number" },
        { title: "R-free", field: "r_free", width: 82, hozAlign: "right", formatter: fixed(3),
          sorter: "number" },
        { title: "Space group", field: "space_group", width: 118, headerFilter: "input",
          formatter: mono },
        { title: "Length", field: "length", width: 82, hozAlign: "right", formatter: mono,
          sorter: "number" },
        { title: "Chains", field: "chains", width: 92, formatter: mono },
        { title: "Organism", field: "organism", width: 168, headerFilter: "input" },
        { title: "Expressed in", field: "host", width: 150, headerFilter: "input", visible: false },
        { title: "UniProt", field: "uniprot", width: 96, formatter: mono, visible: false },
        { title: "Ligands", field: "ligands", width: 168, headerFilter: "input", formatter: mono },
        { title: "Holo", field: "holo", width: 68 },
        { title: "Fusion", field: "fusion", width: 78, visible: false },
        { title: "Deposited", field: "deposited", width: 108, formatter: mono, visible: false },
        { title: "Released", field: "released", width: 104, formatter: mono },
        { title: "Title", field: "title", width: 340, headerFilter: "input", visible: false }
      ]
    });

    // The table is a panel like any other: it binds to the same selection state.
    table.on("rowMouseEnter", function (e, row) { setHot(row.getData().entity_id); });
    table.on("rowMouseLeave", function () { setHot(null); });
    table.on("rowClick", function (e, row) { openCard(row.getData().entity_id); });

    wireExports();
    wireColumnToggle();
  }

  function highlightTableRow(id) {
    if (!table) return;
    try {
      var rows = table.getRows();
      rows.forEach(function (r) {
        r.getElement().classList.toggle("hot", id != null && r.getData().entity_id === id);
      });
    } catch (e) { /* the table is not built yet */ }
  }

  function wireExports() {
    $("expCsv").addEventListener("click", function () {
      table.download("csv", CFG.slug + ".csv");
    });
    $("expJson").addEventListener("click", function () {
      table.download("json", CFG.slug + ".json");
    });
    $("expBib").addEventListener("click", exportBibtex);
  }

  /* One BibTeX record per distinct primary citation in the visible set, not one per entry:
     a family of 1,500 structures is a few hundred papers, and the duplicates are what makes
     a hand-assembled bibliography of a family so tedious. */
  function exportBibtex() {
    var seen = {}, out = [];
    S.members.forEach(function (m) {
      if (!S.visible.has(m.entity_id) || !m.cite_id) return;
      var c = S.family.citations[m.cite_id];
      if (!c) return;
      var key = m.cite_id;
      if (!key || seen[key]) {
        if (seen[key]) seen[key].pdbs.push(m.pdb_id);
        return;
      }
      var rec = { c: c, pdbs: [m.pdb_id] };
      seen[key] = rec;
      out.push(rec);
    });

    var text = out.map(function (rec) {
      var c = rec.c;
      var first = (c.authors && c.authors[0] || "anon").split(",")[0].replace(/\W/g, "");
      var cite = (first + (c.year || "") + (rec.pdbs[0] || "")).toLowerCase();
      var lines = [
        "@article{" + cite + ",",
        "  title   = {" + (c.title || "") + "},",
        "  author  = {" + (c.authors || []).join(" and ") + "},",
        "  journal = {" + (c.journal || "") + "},"
      ];
      if (c.year) lines.push("  year    = {" + c.year + "},");
      if (c.volume) lines.push("  volume  = {" + c.volume + "},");
      if (c.pages) lines.push("  pages   = {" + c.pages + (c.page_last ? "--" + c.page_last : "") + "},");
      if (c.doi) lines.push("  doi     = {" + c.doi + "},");
      lines.push("  note    = {PDB " + rec.pdbs.sort().join(", ") + "}");
      lines.push("}");
      return lines.join("\n");
    }).join("\n\n");

    var header = "% " + (S.family.name || CFG.slug) + " — " + out.length +
      " primary citations across " + S.visible.size + " entities.\n" +
      "% Assembled by CODSWALLOP (codswallop.mdeller.com) from RCSB PDB metadata.\n\n";

    var blob = new Blob([header + text], { type: "application/x-bibtex" });
    var a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = CFG.slug + ".bib";
    a.click();
    URL.revokeObjectURL(a.href);
  }

  function wireColumnToggle() {
    $("colToggle").addEventListener("click", function () {
      var open = document.querySelector(".colmenu");
      if (open) { open.remove(); return; }
      var menu = document.createElement("div");
      menu.className = "colmenu";
      menu.style.cssText = "position:absolute;z-index:20;background:var(--surface);" +
        "border:1px solid var(--line);border-radius:8px;padding:9px 12px;box-shadow:var(--shadow);" +
        "max-height:320px;overflow:auto;font-size:12px";
      table.getColumns().forEach(function (col) {
        var def = col.getDefinition();
        var id = "col-" + def.field;
        var row = document.createElement("label");
        row.style.cssText = "display:flex;gap:7px;align-items:center;padding:3px 0;color:var(--dim);cursor:pointer";
        row.innerHTML = '<input type="checkbox" id="' + id + '"' +
          (col.isVisible() ? " checked" : "") + "> " + esc(def.title);
        row.querySelector("input").addEventListener("change", function () {
          col.toggle();
        });
        menu.appendChild(row);
      });
      var r = this.getBoundingClientRect();
      menu.style.left = (r.left + window.scrollX) + "px";
      menu.style.top = (r.bottom + window.scrollY + 5) + "px";
      document.body.appendChild(menu);
      setTimeout(function () {
        document.addEventListener("click", function once(ev) {
          if (!menu.contains(ev.target)) { menu.remove(); document.removeEventListener("click", once); }
        });
      }, 0);
    });
  }

  /* ==================================================================================
     10. Map / Table segmented control, and the rail
     ================================================================================== */
  function wireView() {
    function set(view) {
      S.view = view;
      $("btnMap").setAttribute("aria-pressed", view === "map" ? "true" : "false");
      $("btnTable").setAttribute("aria-pressed", view === "table" ? "true" : "false");
      $("hero").classList.toggle("table-mode", view === "table");
      $("tableWrap").hidden = view !== "table";
      if (view === "table") {
        buildTable();
        // Tabulator measures itself on construction; if it was built while hidden it comes
        // out zero-height, so nudge it once it is actually on screen.
        requestAnimationFrame(function () { try { table.redraw(true); } catch (e) {} });
      } else if (constellation) {
        requestAnimationFrame(function () { constellation.resize(); constellation.draw(); });
      }
    }
    $("btnMap").addEventListener("click", function () { set("map"); });
    $("btnTable").addEventListener("click", function () { set("table"); });
  }

  function wireRail() {
    var tabs = document.querySelectorAll(".divider");
    tabs.forEach(function (tab) {
      tab.addEventListener("click", function () {
        var id = tab.dataset.section;
        tabs.forEach(function (t) {
          t.setAttribute("aria-selected", t === tab ? "true" : "false");
        });
        document.querySelectorAll(".section").forEach(function (sec) {
          sec.hidden = sec.id !== "sec-" + id;
        });
        if (id === "structures") {
          var pick = $("structurePick");
          // Only once: re-entering the section must not reload 5 MB or reset the camera.
          if (pick && pick.value && !$("structureHost").dataset.loaded) {
            $("structureHost").dataset.loaded = "1";
            loadStructure(pick.value);
          }
        }
        if (id === "overview" && constellation) {
          requestAnimationFrame(function () { constellation.resize(); constellation.draw(); });
        }
        window.scrollTo({ top: 0, behavior: "smooth" });
      });
    });
  }

  /* The "what phase 1 can already tell you" blocks under each placeholder. Real data the
     app has now, shown where the reader is standing when they wonder about it. */
  function renderKnownBlocks() {
    var s = S.family.stats;
    var renderers = {
      ligands: function () {
        if (!s.top_ligands.length) return "<p>No chemical components outside the buffer and ion list.</p>";
        return table2(["Component", "Name", "Entries"], s.top_ligands.map(function (l) {
          return [l.id, l.name || "", commas(l.count)];
        }), [0, 2]);
      },
      crystals: function () {
        if (!s.top_space_groups.length) return "<p>No space groups recorded: nothing here was crystallised.</p>";
        return table2(["Space group", "Entries"], s.top_space_groups.map(function (g) {
          return [g[0], commas(g[1])];
        }), [0, 1]);
      },
      methods: function () {
        return table2(["Method", "Entries"], s.methods.map(function (m) {
          return [m[0], commas(m[1])];
        }), [1]);
      },
      constructs: function () {
        return "<p style='color:var(--dim);font-size:13px;line-height:1.7'>" +
          "The family holds <b style='color:var(--text)'>" + commas(s.constructs) +
          "</b> distinct deposited sequences across " + commas(s.entities) + " entities, of which " +
          "<b style='color:var(--text)'>" + commas(s.fusions) + "</b> carry an apparent fusion " +
          "partner or extension and <b style='color:var(--text)'>" + commas(s.orthologues) +
          "</b> come from an organism other than the seed's. The median construct contains " +
          s.coverage.median_coverage + "% of the seed sequence.</p>";
      },
      domains: function () {
        var f = S.family;
        var rows = (f.pfam || []).map(function (p) { return ["Pfam", p.id, p.name || ""]; })
          .concat((f.interpro || []).map(function (p) { return ["InterPro", p.id, p.name || ""]; }));
        if (!rows.length) return "<p>No Pfam or InterPro annotation on the seed entity.</p>";
        return table2(["Source", "Accession", "Name"], rows, [1]);
      },
      quality: function () {
        var withR = S.members.filter(function (m) { return m.r_free != null && m.r_work != null; });
        if (!withR.length) return "<p>No refinement statistics in this family.</p>";
        var gaps = withR.map(function (m) { return m.r_free - m.r_work; }).sort(function (a, b) { return a - b; });
        var frees = withR.map(function (m) { return m.r_free; }).sort(function (a, b) { return a - b; });
        return table2(["Figure", "Value"], [
          ["Entities with R-free", commas(withR.length)],
          ["Median R-free", frees[Math.floor(frees.length / 2)].toFixed(3)],
          ["Median R-free − R-work", gaps[Math.floor(gaps.length / 2)].toFixed(3)],
          ["Widest R-free − R-work", gaps[gaps.length - 1].toFixed(3)]
        ], [1]);
      }
    };

    document.querySelectorAll("[data-known]").forEach(function (box) {
      var fn = renderers[box.dataset.known];
      if (fn) box.innerHTML = fn();
    });
  }

  function table2(headers, rows, numericCols) {
    return '<table class="minitable"><thead><tr>' +
      headers.map(function (h) { return "<th>" + esc(h) + "</th>"; }).join("") +
      "</tr></thead><tbody>" + rows.slice(0, 20).map(function (r) {
        return "<tr>" + r.map(function (v, i) {
          var cls = numericCols.indexOf(i) >= 0 ? "n" : (i === 0 ? "name" : "");
          return '<td class="' + cls + '">' + esc(v) + "</td>";
        }).join("") + "</tr>";
      }).join("") + "</tbody></table>";
  }

  /* ==================================================================================
     10b. Phase 2: sequences, constructs, domains
     ================================================================================== */

  /* Two curves on one axis. The blue is how many constructs CONTAIN each residue; the green
     is how many actually RESOLVED it. The gap between them is the disorder, and drawing them
     together is the whole point: either curve alone answers a different question, and a
     reader who conflates them designs the wrong construct. */
  function renderCoverage2() {
    var c = S.family.stats.coverage;
    var svg = $("coverage2");
    var depth = c.depth || [];
    if (!depth.length) { svg.innerHTML = ""; return; }
    var seen = c.seen;

    var W = 1000, H = 96, pad = 14, base = H - 16;
    var n = depth.length, maxd = c.max_depth || 1;
    svg.setAttribute("viewBox", "0 0 " + W + " " + H);
    svg.setAttribute("preserveAspectRatio", "none");
    var bw = (W - pad * 2) / n, parts = [];

    for (var i = 0; i < n; i++) {
      var x = (pad + i * bw).toFixed(2), w = Math.max(bw, 0.6).toFixed(2);
      var hD = ((base - 2) * depth[i]) / maxd;
      parts.push('<rect class="depth" x="' + x + '" y="' + (base - hD).toFixed(2) +
        '" width="' + w + '" height="' + Math.max(hD, depth[i] ? 0.6 : 0).toFixed(2) +
        '" opacity="0.42"></rect>');
      if (seen) {
        var hS = ((base - 2) * seen[i]) / maxd;
        parts.push('<rect class="seen" x="' + x + '" y="' + (base - hS).toFixed(2) +
          '" width="' + w + '" height="' + Math.max(hS, 0).toFixed(2) +
          '"><title>residue ' + (i + 1) + ": in " + depth[i] + " constructs, resolved in " +
          seen[i] + "</title></rect>");
      }
    }
    parts.push('<text class="covaxis" x="' + pad + '" y="' + (H - 3) + '">1</text>');
    parts.push('<text class="covaxis" x="' + (W - pad) + '" y="' + (H - 3) +
      '" text-anchor="end">' + n + "</text>");
    svg.innerHTML = parts.join("");

    var legend = '<span><i style="background:var(--accent);opacity:.42"></i>in the construct</span>';
    if (seen) {
      legend += '<span><i style="background:var(--mint)"></i>resolved in the density</span>' +
        "<span>" + c.rarely_resolved_pct + "% of the seed is resolved in under " +
        Math.round(c.rarely_cut * 100) + "% of the constructs that contained it.</span>";
    } else {
      legend += "<span>Density not measured for this family.</span>";
    }
    $("covLegend2").innerHTML = legend;
  }

  function renderDisorderRuns() {
    var c = S.family.stats.coverage;
    var runs = c.disorder_runs || [];
    if (!c.seen) {
      $("disorderRuns").innerHTML = "<p style='color:var(--mute);font-size:12.5px;margin:0'>" +
        "Per-chain density was not available for this family.</p>";
      return;
    }
    if (!runs.length) {
      $("disorderRuns").innerHTML = "<p style='color:var(--mint);font-size:12.5px;margin:0'>" +
        "No region of this protein is systematically unresolved. Every residue that made it " +
        "into a construct was resolved in most of them.</p>";
      return;
    }
    $("disorderRuns").innerHTML = table2(
      ["Residues", "Length", "Resolved in"],
      runs.map(function (r) {
        return [r.start + "\u2013" + r.end, r.length + " aa",
                Math.round(r.resolved_fraction * 100) + "% of the constructs containing them"];
      }), [1]);
  }

  function renderSeedSequence() {
    var f = S.family, c = f.stats.coverage;
    var seq = f.seed_sequence || "";
    if (!seq) { $("seedSeq").innerHTML = ""; return; }
    var d = c.depth || [], se = c.seen;
    var out = [];
    for (var i = 0; i < seq.length; i++) {
      if (i % 60 === 0) {
        out.push('<br><span class="ruler">' + String(i + 1).padStart(5, " ") + " </span>");
      }
      var cls = "ok";
      if (!d[i]) cls = "none";
      else if (se && d[i] && se[i] / d[i] < c.rarely_cut) cls = "poor";
      out.push('<i class="' + cls + '" title="' + (i + 1) + " " + seq[i] +
        (d[i] ? ": in " + d[i] + " constructs" + (se ? ", resolved in " + se[i] : "") : "") +
        '">' + seq[i] + "</i>");
    }
    $("seedSeq").innerHTML = out.join("").replace(/^<br>/, "");
  }

  var constructFilter = "all";

  function renderConstructs() {
    var cs = S.family.constructs || [];
    var f = S.family;
    var shown = cs.filter(function (c) {
      if (constructFilter === "eng") return c.engineered;
      if (constructFilter === "fus") return (c.fusions || []).length;
      return true;
    });
    $("constructSub").textContent = commas(cs.length) + " distinct deposited sequences across " +
      commas(S.members.length) + " entities" +
      (f.stats.constructs_unreferenced ? "; " + f.stats.constructs_unreferenced +
        " have no UniProt reference to diff against" : "");

    if (!shown.length) {
      $("constructTable").innerHTML = '<p class="empty">Nothing matches that filter.</p>';
      return;
    }

    var rows = shown.map(function (c) {
      var badges = [];
      (c.tags || []).forEach(function (t) { badges.push('<span class="badge tag">' + esc(t) + "</span>"); });
      (c.fusions || []).forEach(function (t) {
        badges.push('<span class="badge fus">' + esc(t.split(" (")[0]) + "</span>");
      });
      (c.proteases || []).forEach(function (t) { badges.push('<span class="badge prot">' + esc(t) + "</span>"); });
      if (c.mutation_count) {
        badges.push('<span class="badge mut">' + c.mutation_count + " mut</span>");
      }
      var best = c.best_resolution ? Number(c.best_resolution).toFixed(2) : "\u2014";
      return '<tr data-seq="' + esc(c.seq_id) + '">' +
        '<td class="n">' + commas(c.n_entities) + "</td>" +
        '<td class="n">' + c.length + "</td>" +
        '<td class="n">' + best + "</td>" +
        '<td class="id">' + (c.best_pdb_id
          ? '<a href="https://www.rcsb.org/structure/' + esc(c.best_pdb_id) +
            '" target="_blank" rel="noopener noreferrer">' + esc(c.best_pdb_id) + "</a>"
          : "\u2014") + "</td>" +
        '<td class="n">' + esc(c.uniprot || "\u2014") + "</td>" +
        "<td>" + badges.join("") + '<div class="cs">' + esc(c.summary) + "</div></td>" +
        "</tr>";
    }).join("");

    $("constructTable").innerHTML =
      '<div class="tablewrap"><table class="ctable" id="ctBody"><thead><tr>' +
      "<th>Entities</th><th>Length</th><th>Best (\u00c5)</th><th>Best entry</th>" +
      "<th>Reference</th><th>What was made</th></tr></thead><tbody>" +
      rows + "</tbody></table></div>";
  }

  /* The construct table binds to the same selection state as everything else: hovering a
     construct dims every node on the map that does not use it. This is the rule the design
     hangs on, and a panel that kept its own idea of what is selected would break it. */
  function wireConstructSelection() {
    var box = $("constructTable");
    if (!box) return;
    box.addEventListener("mouseover", function (ev) {
      var tr = ev.target.closest("tr[data-seq]");
      if (tr) setConstruct(tr.getAttribute("data-seq"));
    });
    box.addEventListener("mouseleave", function () { setConstruct(null); });
    box.addEventListener("click", function (ev) {
      var tr = ev.target.closest("tr[data-seq]");
      if (!tr) return;
      // Clicking opens the best entry that used this construct, so the reader lands on a
      // real structure rather than on an abstraction.
      var c = (S.family.constructs || []).find(function (x) {
        return x.seq_id === tr.getAttribute("data-seq");
      });
      var m = c && S.members.find(function (mm) { return mm.pdb_id === c.best_pdb_id; });
      if (m) openCard(m.entity_id);
    });
  }

  function setConstruct(seqId) {
    if (S.construct === seqId) return;
    S.construct = seqId;
    if (!constellation) return;
    if (!seqId) {
      constellation.applyVisible(S.visible);
      return;
    }
    var only = new Set();
    S.members.forEach(function (m) {
      if (m.seq_id === seqId && S.visible.has(m.entity_id)) only.add(m.entity_id);
    });
    constellation.applyVisible(only);
  }

  function wireConstructFilters() {
    var map = { cbAll: "all", cbEng: "eng", cbFus: "fus" };
    Object.keys(map).forEach(function (id) {
      var b = $(id);
      if (!b) return;
      b.addEventListener("click", function () {
        constructFilter = map[id];
        Object.keys(map).forEach(function (o) {
          $(o).setAttribute("aria-pressed", o === id ? "true" : "false");
        });
        renderConstructs();
      });
    });
  }

  function renderDomains() {
    var d = S.family.domains || { sources: [], domains: [] };
    var n = S.family.seed_length || 1;
    $("domainSub").textContent = d.domains.length
      ? d.domains.length + " assignments from " + d.sources.join(", ") +
        ", each seen on at least " + d.support + " chains"
      : "no domain assignments met the support threshold";

    if (!d.domains.length) {
      $("domainRibbon").innerHTML = "<p style='color:var(--mute);font-size:12.5px;margin:0'>" +
        "No domain was assigned consistently enough across this family's chains to draw." +
        "</p>";
      return;
    }
    var bySource = {};
    d.domains.forEach(function (x) { (bySource[x.source] = bySource[x.source] || []).push(x); });

    var html = '<div class="ribbon">';
    Object.keys(bySource).sort().forEach(function (src) {
      html += '<div><div class="tracklabel">' + esc(src) + '</div><div class="track">';
      bySource[src].forEach(function (x) {
        var left = (100 * (x.start - 1) / n), width = (100 * (x.end - x.start + 1) / n);
        html += '<div class="dom ' + esc(src.toLowerCase()) + '" style="left:' +
          left.toFixed(2) + "%;width:" + Math.max(width, 0.5).toFixed(2) + '%" title="' +
          esc((x.name || x.id || "") + " \u00b7 " + x.start + "-" + x.end +
              " \u00b7 " + x.n_chains + " chains") + '">' +
          esc(x.name || x.id || "") + "</div>";
      });
      html += "</div></div>";
    });
    html += '<div class="axis"><span>1</span><span>' + n + "</span></div></div>";
    $("domainRibbon").innerHTML = html;
  }

  function renderOrthologues() {
    var rows = S.family.orthologues || [];
    if (!rows.length) { $("orthoMatrix").innerHTML = ""; return; }
    var body = rows.map(function (o) {
      var res = o.best_resolution ? Number(o.best_resolution).toFixed(2) : "\u2014";
      var cov = o.coverage_pct == null ? 0 : o.coverage_pct;
      return "<tr>" +
        '<td class="org">' + esc(o.organism) + "</td>" +
        '<td class="n">' + commas(o.entries) + "</td>" +
        '<td class="n">' + commas(o.entities) + "</td>" +
        '<td class="n">' + res + "</td>" +
        '<td class="n">' + commas(o.holo) + "</td>" +
        '<td><div class="bartrack"><i class="bar" style="width:' + cov.toFixed(0) +
          '%"></i></div></td>' +
        '<td class="n">' + cov + "%</td>" +
        '<td class="n">' + esc((o.accessions || []).join(" ")) + "</td>" +
        "</tr>";
    }).join("");
    $("orthoMatrix").innerHTML = '<table class="omatrix"><thead><tr>' +
      "<th>Organism</th><th>Entries</th><th>Entities</th><th>Best \u00c5</th><th>Holo</th>" +
      "<th>Seed coverage</th><th></th><th>Accessions</th></tr></thead><tbody>" +
      body + "</tbody></table>";
  }

  /* ==================================================================================
     10c. Phase 3 (metadata half): ligands, crystals, quality
     ================================================================================== */
  var ligandFilter = "real";
  var DESIGNED = { ligand: 1, cofactor: 1 };

  function renderLigands() {
    var L = S.family.ligands || { components: [], by_class: {} };
    var comps = L.components.filter(function (c) {
      return ligandFilter === "all" || DESIGNED[c.klass];
    });
    var counts = Object.keys(L.by_class || {}).sort(function (a, b) {
      return L.by_class[b] - L.by_class[a];
    }).map(function (k) { return k + " " + L.by_class[k]; }).join(" · ");
    $("ligandSub").textContent = commas(L.n) + " distinct components · " + counts;

    if (!comps.length) {
      $("ligandTable").innerHTML = '<p class="empty">No components in that class.</p>';
      return;
    }
    var rows = comps.slice(0, 400).map(function (c) {
      return "<tr>" +
        '<td class="n"><a href="https://www.rcsb.org/ligand/' + esc(c.id) +
          '" target="_blank" rel="noopener noreferrer">' + esc(c.id) + "</a></td>" +
        '<td class="n">' + commas(c.count) + "</td>" +
        '<td class="n">' + (c.best_resolution ? Number(c.best_resolution).toFixed(2) : "—") + "</td>" +
        '<td><span class="klass ' + esc(c.klass.split("/")[0]) + '">' + esc(c.klass) + "</span></td>" +
        "<td>" + esc(c.name || "") +
          (c.smiles ? '<span class="smiles">' + esc(c.smiles) + "</span>" : "") + "</td>" +
        '<td class="n">' + esc(c.formula || "") + "</td>" +
        "</tr>";
    }).join("");
    $("ligandTable").innerHTML =
      '<div class="tablewrap"><table class="ctable"><thead><tr>' +
      "<th>CCD</th><th>Entries</th><th>Best Å</th><th>Class</th><th>Name</th><th>Formula</th>" +
      "</tr></thead><tbody>" + rows + "</tbody></table></div>" +
      '<p class="caveat">Metals are classified as ions. Whether a metal is structural or ' +
      "catalytic depends on the protein, not the component: the zinc in carbonic anhydrase " +
      "is a catalytic cofactor filed here as an ion. " +
      "&ldquo;Ligand-bound&rdquo; on the overview counts ligands and cofactors only.</p>";
  }

  function renderCrystals() {
    var c = S.family.crystals || { n: 0 };
    if (!c.n) {
      $("crystalSub").textContent = "no crystallisation conditions recorded";
      $("crystalTables").innerHTML = '<p class="empty">Nothing here was crystallised.</p>';
      return;
    }
    $("crystalSub").textContent = commas(c.n) + " entries carry a condition, " +
      commas(c.n_parsed) + " parsed · pH " + c.ph_min + "–" + c.ph_max +
      " (median " + c.ph_median + ")" +
      (c.temp_median != null ? " · median " + c.temp_median + "°C" : "");

    // pH against precipitant, coloured by resolution. The question is "what pH and what
    // precipitant produced the good crystals", so resolution has to be the colour.
    var pts = c.points || [];
    var svg = $("crystalScatter");
    if (!pts.length) { svg.innerHTML = ""; }
    else {
      var W = 1000, H = 280, padL = 46, padB = 34, padT = 12, padR = 12;
      var cats = [];
      pts.forEach(function (p) { if (cats.indexOf(p.precipitant) < 0) cats.push(p.precipitant); });
      cats = cats.slice(0, 14);
      var res = pts.map(function (p) { return p.resolution; }).filter(function (r) { return r != null; });
      var rMin = Math.min.apply(null, res), rMax = Math.max.apply(null, res);
      var phLo = Math.floor(Math.min.apply(null, pts.map(function (p) { return p.ph; })));
      var phHi = Math.ceil(Math.max.apply(null, pts.map(function (p) { return p.ph; })));
      var X = function (ph) { return padL + (W - padL - padR) * (ph - phLo) / Math.max(1, phHi - phLo); };
      var Y = function (cat) {
        var i = cats.indexOf(cat); if (i < 0) i = cats.length;
        return padT + (H - padT - padB) * (i + 0.5) / (cats.length + 1);
      };
      var out = [];
      for (var t = phLo; t <= phHi; t++) {
        out.push('<line x1="' + X(t) + '" y1="' + padT + '" x2="' + X(t) + '" y2="' + (H - padB) +
          '" stroke="var(--line)" stroke-width="0.6"/>');
        out.push('<text class="covaxis" x="' + X(t) + '" y="' + (H - padB + 14) +
          '" text-anchor="middle">pH ' + t + "</text>");
      }
      cats.forEach(function (cat) {
        out.push('<text class="covaxis" x="4" y="' + (Y(cat) + 3) + '">' +
          esc(cat.length > 14 ? cat.slice(0, 13) + "…" : cat) + "</text>");
      });
      pts.forEach(function (p) {
        if (p.ph == null) return;
        // Better resolution is a hotter dot: mint at the good end, oxide at the poor end.
        var f = (p.resolution == null || rMax === rMin) ? 0.5
              : (p.resolution - rMin) / (rMax - rMin);
        var col = f < 0.34 ? "var(--mint)" : (f < 0.67 ? "var(--cyan)" : "var(--oxide)");
        out.push('<circle class="scatterdot" cx="' + X(p.ph).toFixed(1) + '" cy="' +
          (Y(p.precipitant) + (Math.random() - 0.5) * 9).toFixed(1) +
          '" r="3.1" fill="' + col + '" fill-opacity="0.75"><title>' + esc(p.pdb_id) +
          " · pH " + p.ph + (p.resolution ? " · " + p.resolution + " Å" : "") +
          " · " + esc(p.precipitant) + "</title></circle>");
      });
      svg.setAttribute("viewBox", "0 0 " + W + " " + H);
      svg.setAttribute("preserveAspectRatio", "none");
      svg.innerHTML = out.join("");
      $("crystalLegend").innerHTML =
        '<span><i style="background:var(--mint)"></i>best resolution</span>' +
        '<span><i style="background:var(--cyan)"></i>middle</span>' +
        '<span><i style="background:var(--oxide)"></i>poorest</span>' +
        "<span>" + pts.length + " entries with a recorded pH. Vertical jitter only.</span>";
    }

    function condTable(title, rows) {
      if (!rows || !rows.length) return "";
      return "<div><h3 style='font-size:10px;text-transform:uppercase;letter-spacing:.12em;" +
        "color:var(--mute);margin:0 0 8px'>" + esc(title) + "</h3>" +
        table2(["Component", "Entries", "Best Å"], rows.map(function (r) {
          return [r.name, commas(r.count),
                  r.best_resolution ? Number(r.best_resolution).toFixed(2) : "—"];
        }), [1, 2]) + "</div>";
    }
    $("crystalTables").innerHTML = '<div class="condgrid">' +
      condTable("Precipitants", c.precipitants) +
      condTable("Buffers", c.buffers) +
      condTable("Methods", c.methods) +
      condTable("Additives", c.additives) + "</div>";
  }

  function exportCrystals() {
    var c = S.family.crystals || {};
    var rows = [["pdb_id", "ph", "temp_c", "precipitant", "buffer", "resolution"]];
    (c.points || []).forEach(function (p) {
      rows.push([p.pdb_id, p.ph, p.temp_c == null ? "" : p.temp_c,
                 p.precipitant, p.buffer || "", p.resolution == null ? "" : p.resolution]);
    });
    var csv = rows.map(function (r) {
      return r.map(function (v) {
        var s = String(v == null ? "" : v);
        return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
      }).join(",");
    }).join("\n");
    var a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([csv], { type: "text/csv" }));
    a.download = CFG.slug + "-conditions.csv";
    a.click();
    URL.revokeObjectURL(a.href);
  }

  function renderQuality() {
    var q = S.family.quality || { n: 0, rows: [] };
    if (!q.n) {
      $("qualitySub").textContent = "no validation reports for this family";
      $("qualityTable").innerHTML = '<p class="empty">Nothing to triage.</p>';
      return;
    }
    $("qualitySub").innerHTML = commas(q.n) + " validated · " +
      '<span class="verdict ok"></span>' + commas(q.ok) + " clean · " +
      '<span class="verdict check"></span>' + commas(q.check) + " worth a look · " +
      '<span class="verdict poor"></span>' + commas(q.poor) + " flagged twice or more · " +
      commas(q.with_sf) + " with structure factors";

    var rows = q.rows.map(function (r) {
      var n = function (v, dp) {
        return v == null ? "—" : Number(v).toFixed(dp === undefined ? 2 : dp);
      };
      return "<tr>" +
        '<td class="n"><span class="verdict ' + r.verdict + '"></span>' +
          '<a href="https://www.rcsb.org/structure/' + esc(r.pdb_id) +
          '" target="_blank" rel="noopener noreferrer">' + esc(r.pdb_id) + "</a></td>" +
        '<td class="n">' + n(r.resolution) + "</td>" +
        '<td class="n">' + n(r.clashscore, 1) + "</td>" +
        '<td class="n">' + n(r.rsrz, 1) + "</td>" +
        '<td class="n">' + n(r.rama, 1) + "</td>" +
        '<td class="n">' + n(r.rota, 1) + "</td>" +
        '<td class="n">' + n(r.r_gap, 3) + "</td>" +
        '<td class="n">' + n(r.completeness, 1) + "</td>" +
        '<td class="n">' + (r.has_sf ? "yes" : "no") + "</td>" +
        "<td>" + esc((r.flags || []).join(", ")) + "</td>" +
        "</tr>";
    }).join("");
    $("qualityTable").innerHTML =
      '<div class="tablewrap"><table class="ctable"><thead><tr>' +
      "<th>Entry</th><th>Res Å</th><th>Clash</th><th>RSRZ %</th><th>Rama %</th>" +
      "<th>Rota %</th><th>R gap</th><th>Compl %</th><th>SF</th><th>Flagged</th>" +
      "</tr></thead><tbody>" + rows + "</tbody></table></div>" +
      '<p class="caveat">Worst first. Thresholds are the wwPDB report&rsquo;s own: ' +
      "clashscore over 20, RSRZ or rotamer outliers over 5&nbsp;%, Ramachandran outliers " +
      "over 0.5&nbsp;%, R-free minus R-work over 0.07. Family medians: clashscore " +
      (q.median_clashscore == null ? "—" : q.median_clashscore) + ", RSRZ " +
      (q.median_rsrz == null ? "—" : q.median_rsrz) + "&nbsp;%, R gap " +
      (q.median_r_gap == null ? "—" : q.median_r_gap) + ".</p>";
  }

  function wirePhase3() {
    var lg = { lgAll: "all", lgReal: "real" };
    Object.keys(lg).forEach(function (id) {
      var b = $(id);
      if (!b) return;
      b.addEventListener("click", function () {
        ligandFilter = lg[id];
        Object.keys(lg).forEach(function (o) {
          $(o).setAttribute("aria-pressed", o === id ? "true" : "false");
        });
        renderLigands();
      });
    });
    var ex = $("expCrystals");
    if (ex) ex.addEventListener("click", exportCrystals);
    var cutter = $("tmCut");
    if (cutter) cutter.addEventListener("input", renderHeatmap);
    var sp = $("btnSuperpose");
    if (sp) sp.addEventListener("click", superposeFamily);
    var af = $("btnAlphaFold");
    if (af) af.addEventListener("click", addAlphaFold);
  }

  /* ---- conservation, the logo, and the engineered positions ----------------------- */
  function renderConservation() {
    var m = S.family.msa;
    if (!m || m.too_long || !m.columns || !m.columns.length) {
      $("consSub").textContent = m && m.too_long
        ? "seed too long for a per-residue view" : "not computed for this family";
      $("conservation").innerHTML = "";
      $("engineered").innerHTML = "";
      return;
    }
    $("consSub").textContent = m.n_sequences + " distinct constructs aligned to the seed · " +
      "mean conservation " + m.mean_conservation + " · " + m.conserved.length +
      " positions above " + 0.97;

    var W = 1000, H = 96, pad = 14, base = H - 16;
    var cols = m.columns, n = cols.length;
    var bw = (W - pad * 2) / n, parts = [];
    for (var i = 0; i < n; i++) {
      var c = cols[i];
      if (c.conservation == null) continue;
      var h = (base - 2) * c.conservation;
      // Above the threshold the bar goes mint: those are the positions worth reading.
      var cls = c.conservation >= 0.97 ? "cons hi" : "cons";
      parts.push('<rect class="' + cls + '" x="' + (pad + i * bw).toFixed(2) + '" y="' +
        (base - h).toFixed(2) + '" width="' + Math.max(bw, 0.6).toFixed(2) + '" height="' +
        Math.max(h, 0).toFixed(2) + '"><title>' + c.seed + (i + 1) + ": conservation " +
        c.conservation + ", " + (c.top[0] ? c.top[0].aa + " " +
        Math.round(c.top[0].f * 100) + "%" : "") + "</title></rect>");
    }
    parts.push('<text class="covaxis" x="' + pad + '" y="' + (H - 3) + '">1</text>');
    parts.push('<text class="covaxis" x="' + (W - pad) + '" y="' + (H - 3) +
      '" text-anchor="end">' + n + "</text>");
    var svg = $("conservation");
    svg.setAttribute("viewBox", "0 0 " + W + " " + H);
    svg.setAttribute("preserveAspectRatio", "none");
    svg.innerHTML = parts.join("");
    $("consLegend").innerHTML =
      '<span><i style="background:var(--mint)"></i>conserved (&ge;0.97)</span>' +
      '<span><i style="background:var(--violet)"></i>variable</span>' +
      "<span>Normalised Shannon entropy over the 20 residues, weighted by how many " +
      "entities used each construct. This is what was <em>deposited</em>, not what evolved: " +
      "a heavily engineered position reads as variable.</span>";

    renderEngineered(m);
  }

  /* A sequence logo for one column: residues stacked by frequency, tallest first. */
  function logoFor(col) {
    if (!col || !col.top || !col.top.length) return "";
    return '<span class="logo">' + col.top.map(function (t) {
      var px = Math.max(6, Math.round(t.f * 30));
      return '<i class="aa-' + esc(t.group) + '" style="font-size:' + px + "px;line-height:" +
        px + 'px" title="' + esc(t.aa) + " " + Math.round(t.f * 100) + '%">' +
        esc(t.aa) + "</i>";
    }).join("") + "</span>";
  }

  function renderEngineered(m) {
    var eng = m.engineered || [];
    if (!eng.length) {
      $("engineered").innerHTML = "<p style='color:var(--mute);font-size:12.5px;margin:0'>" +
        "No position in this family carries a substantial minority substitution.</p>";
      return;
    }
    var byPos = {};
    m.columns.forEach(function (c) { byPos[c.pos] = c; });
    var rows = eng.slice(0, 150).map(function (e) {
      return "<tr>" +
        '<td class="n">' + esc(e.seed) + e.pos + "</td>" +
        '<td class="n">' + Math.round(e.substituted * 1000) / 10 + "%</td>" +
        '<td class="n">' + (e.conservation == null ? "—" : e.conservation) + "</td>" +
        "<td>" + e.variants.map(function (v) {
          return '<span class="badge mut">' + esc(e.seed) + e.pos + esc(v.aa) + " " +
            Math.round(v.f * 1000) / 10 + "%</span>";
        }).join("") + "</td>" +
        "<td>" + logoFor(byPos[e.pos]) + "</td>" +
        "</tr>";
    }).join("");
    $("engineered").innerHTML =
      '<div class="tablewrap"><table class="ctable"><thead><tr>' +
      "<th>Position</th><th>Substituted</th><th>Conservation</th><th>Variants</th>" +
      "<th>Logo</th></tr></thead><tbody>" + rows + "</tbody></table></div>";
  }

  /* The map says which of the two things it is showing. A placeholder that looks like a
     measurement is worse than no map, and an embedding that goes unlabelled wastes the
     one thing that makes it worth the compute. */
  function renderFieldNote(map) {
    var el = $("fieldNote");
    if (!el) return;
    if (map && map.embedded) {
      el.innerHTML = "<b>Structural embedding.</b> Position is multidimensional scaling of " +
        "the pairwise TM-score matrix over " + map.n_representatives +
        " representative structures (" + commas(map.n_pairs) + " alignments, median TM " +
        map.median_tm + "), so distance is structural distance and the 0.5 contour is the " +
        "conventional same-fold boundary. Size is 1/resolution; colour is method; an amber " +
        "halo means ligand-bound." +
        (map.approximated
          ? " " + commas(map.approximated) + " of " + commas(map.nodes.length) +
            " entries use a construct that was not among the representatives and are placed " +
            "at the nearest one by identity."
          : "");
    } else {
      el.innerHTML = "<b>Placeholder layout.</b> Outward is decreasing identity to the seed, " +
        "ranked; sector is source organism; size is 1/resolution; colour is method; an amber " +
        "halo means ligand-bound. This family has no structural embedding yet: run " +
        "<code>CODSWALLOP.py embed</code> on a workstation and the positions become a real " +
        "TM-score embedding.";
    }
  }

  /* ---- the TM-score matrix --------------------------------------------------------- */
  var tmState = null;

  function renderHeatmap() {
    var m = S.family.map;
    var host = $("tmHeatmap");
    if (!host) return;
    if (!m || !m.embedded || !m.tm || !m.tm.length) {
      $("tmSub").textContent = "no structural embedding for this family yet";
      host.innerHTML = '<p class="empty">Run <code>CODSWALLOP.py embed</code> on a ' +
        "workstation and the matrix appears here.</p>";
      $("tmClusters").innerHTML = "";
      return;
    }
    var cut = Number($("tmCut").value) / 100;
    $("tmCutVal").textContent = cut.toFixed(2);
    $("tmSub").textContent = m.n_representatives + " representatives · " +
      commas(m.n_pairs) + " alignments · median TM " + m.median_tm;

    tmState = TmHeatmap.draw(host, m.tm, m.representatives, {
      cutHeight: cut,
      onPickCluster: function (reps) { selectCluster(reps); },
    });
    renderClusterList(m, cut);
    $("tmLegend").innerHTML =
      '<span><i style="background:var(--mint)"></i>TM 1.0, same structure</span>' +
      '<span><i style="background:var(--accent)"></i>TM 0.6</span>' +
      '<span><i style="background:var(--sky);border:1px solid var(--line)"></i>TM &le;0.2</span>' +
      "<span>Average linkage, cut at 1&nbsp;&minus;&nbsp;TM = " + cut.toFixed(2) +
      ". TM 0.5 is the conventional same-fold boundary and is the wrong place to cut here: " +
      "every member of a family is the same fold, so a cut there returns one cluster for " +
      "almost any family. The informative range is much tighter, where conformational " +
      "states and fusion constructs separate. Drag it to see whether a grouping is real or " +
      "an artefact of where the line sits.</span>";
  }

  function renderClusterList(m, cut) {
    if (!tmState) return;
    var box = $("tmClusters");
    box.innerHTML = "<h4>" + tmState.groups.length + " clusters at this cut</h4>" +
      tmState.groups.map(function (g, gi) {
        var ids = g.map(function (l) { return m.representatives[l].pdb_id; });
        var entities = g.reduce(function (a, l) {
          return a + (m.representatives[l].n_entities || 0);
        }, 0);
        return '<div class="tmgroup" data-cluster="' + gi + '">' +
          '<div class="n">' + g.length + " structure" + (g.length === 1 ? "" : "s") +
          " · " + commas(entities) + " entities</div>" +
          '<div class="ids">' + esc(ids.slice(0, 12).join(" ")) +
          (ids.length > 12 ? " +" + (ids.length - 12) : "") + "</div></div>";
      }).join("");

    box.querySelectorAll(".tmgroup").forEach(function (el) {
      el.addEventListener("click", function () {
        var gi = Number(el.dataset.cluster);
        var already = el.classList.contains("on");
        box.querySelectorAll(".tmgroup").forEach(function (o) { o.classList.remove("on"); });
        if (already) { selectCluster(null); return; }
        el.classList.add("on");
        selectCluster(tmState.groups[gi].map(function (l) { return m.representatives[l]; }));
      });
    });
  }

  /* A cluster the matrix found becomes a filter, like any other: the map, the entry list
     and the table all narrow to the entities whose construct is in it. This is what "the
     clusters become first-class filters" has to mean if the shared selection state is
     load-bearing. */
  function selectCluster(reps) {
    S.cluster = reps ? new Set(reps.map(function (r) { return r.seq_id; })) : null;
    apply();
  }

  /* ---- the Structures panel ------------------------------------------------------- */
  function renderStructures() {
    var pick = $("structurePick");
    if (!pick) return;
    // Best resolution first: the entry someone should look at before any other.
    var shown = S.members.filter(function (m) { return S.visible.has(m.entity_id); })
                         .slice(0, 300);
    $("structureSub").textContent = shown.length
      ? "best-resolution entries first" : "nothing passes the current filters";
    pick.innerHTML = shown.map(function (m) {
      return '<option value="' + esc(m.pdb_id) + '">' + esc(m.pdb_id) +
        (m.resolution ? " · " + Number(m.resolution).toFixed(2) + " Å" : " · " + esc(m.method)) +
        " · " + esc((m.description || "").slice(0, 40)) + "</option>";
    }).join("");
    if (!pick.dataset.wired) {
      pick.dataset.wired = "1";
      pick.addEventListener("change", function () { loadStructure(pick.value); });
    }
  }

  var currentViewer = null;

  function loadStructure(pdbId) {
    if (!pdbId) return;
    window.CodswallopViewer.show($("structureHost"), pdbId).then(function (v) {
      currentViewer = v;
      $("superposeNote").textContent = "";
    });
  }

  /* Superpose the representatives of the current TM cluster, or the family's top ones.
     Capped hard: each structure is a separate download and a separate Mol* state tree, and
     a dozen is already a slow load on a laptop. */
  // Six, not ten. Each structure is a separate download, parse and Mol* state-tree
  // commit, done sequentially because the builder cannot be raced; ten took over a
  // minute on a warm connection, which reads as a hang.
  var SUPERPOSE_MAX = 6;

  function superposeFamily() {
    var m = S.family.map;
    if (!m || !m.embedded) {
      $("superposeNote").textContent =
        "Superposition needs the structural embedding: run CODSWALLOP.py embed first.";
      return;
    }
    // If a cluster is selected, superpose that; otherwise the most-used constructs.
    var reps = m.representatives.filter(function (r) {
      return !S.cluster || S.cluster.has(r.seq_id);
    });
    reps = reps.filter(function (r) { return r.transform; })
               .sort(function (a, b) { return b.n_entities - a.n_entities; })
               .slice(0, SUPERPOSE_MAX);
    // The reference must be among the structures actually loaded: every transform maps onto
    // its frame, so superposing without it anchors the pile on a structure that is not there.
    if (m.reference && !reps.some(function (r) { return r.pdb_id === m.reference; })) {
      var ref = m.representatives.find(function (r) { return r.pdb_id === m.reference; });
      if (ref && ref.transform) reps = [ref].concat(reps.slice(0, SUPERPOSE_MAX - 1));
    }
    if (reps.length < 2) {
      $("superposeNote").textContent = "Fewer than two alignable structures here.";
      return;
    }
    // Cool colours only: these are data, and the brass is reserved for the drawer.
    var palette = ["--cyan", "--mint", "--accent", "--violet", "--amber"];
    var entries = reps.map(function (r, i) {
      return {
        pdb_id: r.pdb_id, transform: r.transform,
        colour: cssColour(palette[i % palette.length]),
      };
    });
    $("superposeNote").textContent = "Superposing " + reps.length + " structures onto " +
      (m.reference || reps[0].pdb_id) + "…";
    window.CodswallopViewer.superpose($("structureHost"), entries).then(function (v) {
      currentViewer = v;
      $("superposeNote").textContent = reps.length + " structures on " +
        (m.reference || reps[0].pdb_id) + " (" +
        reps.map(function (r) { return r.pdb_id; }).join(" ") + ")";
    });
  }

  /** A CSS custom property as the {r,g,b} Mol* uniform colouring wants. */
  function cssColour(varName) {
    var hex = getComputedStyle(document.documentElement).getPropertyValue(varName).trim();
    var m = /^#?([0-9a-f]{6})$/i.exec(hex);
    return m ? parseInt(m[1], 16) : 0x5b8cff;
  }

  function addAlphaFold() {
    var acc = S.family.seed && /^[A-Z0-9]{6,10}$/i.test(S.family.seed) ? S.family.seed : null;
    if (!acc) {
      // Fall back to the accession most members carry, since a family seeded from a PDB ID
      // has no accession of its own.
      var counts = {};
      S.members.forEach(function (m) {
        if (m.uniprot) counts[m.uniprot] = (counts[m.uniprot] || 0) + 1;
      });
      acc = Object.keys(counts).sort(function (a, b) { return counts[b] - counts[a]; })[0];
    }
    if (!acc) {
      $("superposeNote").textContent = "No UniProt accession for this family.";
      return;
    }
    if (!currentViewer) {
      $("superposeNote").textContent = "Load a structure first.";
      return;
    }
    $("superposeNote").textContent = "Fetching the AlphaFold model for " + acc + "…";
    window.CodswallopViewer.addAlphaFold(currentViewer, acc).then(function () {
      $("superposeNote").textContent = "AlphaFold model for " + acc +
        " added, coloured by pLDDT. It sits in its own frame: it is not superposed.";
    }).catch(function (err) {
      $("superposeNote").textContent = "No AlphaFold model for " + acc + " (" + err.message + ").";
    });
  }

  /* ==================================================================================
     11. Boot
     ================================================================================== */
  function paint(fam) {
    S.family = fam;
    S.members = fam.members || [];

    $("filing").hidden = true;

    renderHeader();
    renderStats();
    renderFilters();
    renderCoverage();
    renderCoverage2();
    renderDisorderRuns();
    renderSeedSequence();
    renderConservation();
    renderHeatmap();
    renderConstructs();
    renderDomains();
    renderOrthologues();
    renderLigands();
    renderCrystals();
    renderQuality();
    renderKnownBlocks();

    recompute();

    constellation = new Constellation($("constellation"), {
      onHover: setHot,
      onPick: openCard
    });
    renderFieldNote(fam.map);
    constellation.render(fam.map, S.members);

    subscribe(function () {
      renderList();
      // Here rather than in the pre-render block above: S.visible is empty until
      // recompute() runs, so rendering the picker earlier produced an empty <select> and
      // the Structures panel silently never loaded anything. It also has to follow the
      // filters, which is the same reason.
      renderStructures();
      if (constellation) constellation.applyVisible(S.visible);
      if (table) table.replaceData(tableData());
      $("filterSum").innerHTML = "<b>" + commas(S.visible.size) + "</b> of " +
        commas(S.members.length) + " entities";
      $("heroSub").textContent = fam.stats.entries.toLocaleString("en-GB") +
        " entries · " + fam.stats.constructs + " distinct constructs · " +
        (fam.stats.identity_min != null ? fam.stats.identity_min + "–" +
          fam.stats.identity_max + "% identity to seed" : "");
    });

    notify("boot");
  }

  function fail(message) {
    $("filing").innerHTML = '<span class="plate">Empty</span>' +
      '<p>' + esc(message) + "</p>" +
      '<p style="margin-top:16px"><a href="/">← Ask for something else</a></p>';
    $("filing").hidden = false;
  }

  function boot() {
    wireFilters();
    wireList();
    wireView();
    wireRail();
    wireConstructFilters();
    wireConstructSelection();
    wirePhase3();

    $("cardClose").addEventListener("click", closeCard);
    $("scrim").addEventListener("click", closeCard);
    document.addEventListener("keydown", function (ev) {
      if (ev.key === "Escape" && S.picked) closeCard();
    });

    if (CFG.family) { paint(CFG.family); return; }

    // Cold: the shell is already on screen in its filing state, so this is the only wait.
    var steps = ["resolving the query…", "searching the PDB by sequence…",
      "fetching entry metadata…", "laying out the family…"];
    var i = 0;
    var tick = setInterval(function () {
      i = Math.min(i + 1, steps.length - 1);
      $("filingStep").textContent = steps[i];
    }, 2600);
    $("filingStep").textContent = steps[0];

    var url = CFG.apiUrl + (CFG.query ? "?q=" + encodeURIComponent(CFG.query) : "");
    fetch(url, { headers: { "Accept": "application/json" } })
      .then(function (r) { return r.json().then(function (b) { return { ok: r.ok, body: b }; }); })
      .then(function (res) {
        clearInterval(tick);
        if (!res.ok || res.body.error) {
          fail(res.body.error || "The family could not be assembled.");
          return;
        }
        paint(res.body);
      })
      .catch(function (err) {
        clearInterval(tick);
        fail("The request failed: " + err.message);
      });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", boot);
  } else {
    boot();
  }
})();
