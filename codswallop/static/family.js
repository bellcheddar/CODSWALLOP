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
      {
        k: "Thinly covered", flag: true,
        v: c.thin_pct + '%<small> of seed</small>',
        title: "Seed residues present in fewer than " + c.thin_cut + " of the family's " +
          "constructs (" + c.thin + " of " + c.length + " residues). Phase 2 adds the " +
          "harder question: which of them has nobody ever seen density for."
      },
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
    renderKnownBlocks();

    recompute();

    constellation = new Constellation($("constellation"), {
      onHover: setHot,
      onPick: openCard
    });
    constellation.render(fam.map, S.members);

    subscribe(function () {
      renderList();
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
