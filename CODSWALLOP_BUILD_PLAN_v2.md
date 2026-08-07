# CODSWALLOP — Build Plan v2

**C**omparative **O**verview of **D**omains, **S**tructures, **W**hole-families,
**A**lignments, **L**igands, **L**ineages, **O**rganisms & **P**DB entries

> A structural biologist has just been handed a protein family and told to "get up to
> speed". CODSWALLOP is the tool that does the fortnight of manual PDB trawling in
> ninety seconds: every entry, every construct, every ligand, every crystallisation
> condition, cross-referenced and interactive.

- **Repo:** `github.com/bellcheddar/CODSWALLOP`
- **Live at:** `codswallop.mdeller.com`
- **Stack:** Python 3.11 + Flask + gunicorn + nginx + SQLite, vanilla JS front end
  (Mol\*, Plotly.js, D3, Tabulator). Same droplet pattern as AlphaFraud.
- **Author:** Marc C. Deller, D.Phil. (marc@marcdeller.com)
- **Reference prototype:** `prototypes/codswallop_unified_theme_v1.html`

**Changed in v2:** the three competing design concepts have been cut. The Drawer +
Constellation hybrid is now the single agreed direction, specified in full in section 1.

---

## 0. Why this is not just another PDB mirror

The RCSB and PDBe landing pages are *entry-centric*. They answer "what is 4XYZ?".
Nobody answers the question a working structural biologist actually asks:

> "What do the 87 structures of this family collectively tell me, and what should
> I do differently?"

CODSWALLOP is **family-centric**. Its five differentiators:

| # | Differentiator | Why it matters |
|---|---|---|
| 1 | **Construct diffing** | Aligns every deposited SEQRES against the UniProt canonical to expose tags, fusion partners, truncations, surface-entropy-reduction and thermostabilising mutations. This is what actually gets a crystal, and it is buried in mmCIF nobody reads. |
| 2 | **Global cross-filtering** | One selection state shared by every panel. Click a ligand and the map, the entry list, the crystallisation scatter and the coverage plot all filter together. |
| 3 | **Crystallisation intelligence** | Parses `_exptl_crystal_grow` across the whole family into structured conditions, so you can see what actually worked. Feeds directly into the Top96 crystallisation-predictor project. |
| 4 | **Disorder / coverage census** | Which residues has *nobody* ever seen? Family-wide missing-density map answers a construct-design question in one glance. |
| 5 | **Interaction fingerprints** | PLIP run family-wide, so binding-site contacts are comparable across entries rather than one-off. |

---

## 1. Design direction: the Drawer, with the Constellation as its hero

**Agreed direction.** An archival drawer shell provides navigation and identity; a
similarity map provides the hero view inside it. Neither works as well alone: the drawer
alone is a nicely dressed table, and the map alone is a bad first five seconds for
anyone who just wants the resolutions.

### 1.1 Shell

- **Landing:** a closed drawer front, brass plate, one search field. Submit and the
  drawer slides open (~680 ms, `cubic-bezier(.62,.02,.26,1)`) to reveal the family.
- **Navigation:** a persistent vertical rail of **tabbed dividers**, not a top nav.
  Ten sections fit comfortably where a horizontal bar would wrap:
  Overview · Sequences · Constructs · Domains · Structures · Ligands · Contacts ·
  Crystals · Quality · Topology.
- **The pip on each divider is the build phase**, not decoration. Sections that do not
  exist yet open a placeholder naming their phase and what will land there. The rail
  doubles as the roadmap, and the navigation never needs rebuilding as phases ship.
- **Detail view:** clicking any entry, from anywhere, flips out an **index card**
  (slide-over from the right, with scrim, Escape to dismiss). Never a page navigation.
  The family context is never lost.
- **Family header:** name, Pfam/fold/organism metadata, and a rubber-stamped `FILED`
  mark carrying the entry count. Stat strip beneath it, with the "residues never seen"
  figure deliberately flagged in oxide red.

### 1.2 Hero: the family map

- D3-style radial layout of every entry. **Position** = structural similarity (a real
  embedding from the pairwise TM-score matrix once Phase 3 lands, a placeholder cluster
  layout before that). **Radius** = 1/resolution. **Fill** = method. **Amber halo** =
  ligand-bound. **Edges** between entries above the identity threshold.
- Cluster labels are printed on the field, so conformational states (apo open, holo
  closed, thermophilic orthologues) are legible before you read a single number.
- An **entry list sits permanently beside the map**, so numbers are always one glance
  away. A Map/Table segmented control in the panel header swaps to a full sortable table
  for anyone who wants nothing else.

### 1.3 Signature: cross-highlight

Hovering or focusing any row in the entry list pulses that entry's node brass in the
map, and hovering a node highlights its row. One state, two representations. **This is
the single element the interface is remembered by, and it is what makes it read as one
instrument rather than two panels sharing a page.** It also generalises: every panel
added in later phases binds to the same selection state.

### 1.4 Theme: lamp-lit archive, dark by default

Deep blue-slate, never black. The archival identity survives the dark treatment by
moving to brass: anything physical (drawer front, dividers, plate, pins, stamps) is warm;
all data is cool. A light variant is available via a toggle and honours
`prefers-color-scheme`.

| Token | Dark | Role |
|---|---|---|
| `--bg` | `#131924` | Room |
| `--bg-2` | `#182031` | Drawer interior, rail |
| `--sky` | `#0F1522` | Constellation field |
| `--surface` | `#1F2938` | Cards, panels |
| `--surface-2` | `#273346` | Raised, panel headers |
| `--line` | `#2E3A4E` | Borders |
| `--text` / `--dim` / `--mute` | `#E4EAF6` / `#93A2BC` / `#65748E` | Type hierarchy |
| `--brass` / `--brass-dk` | `#C9A063` / `#8A6D40` | Anything physical |
| `--accent` | `#5B8CFF` | Brand blue, lifted for dark contrast |
| `--cyan` / `--amber` / `--violet` | `#4FD1E3` / `#F5B93F` / `#9B7BE8` | X-ray / cryo-EM / NMR |
| `--oxide` | `#E0685C` | FILED stamp, warnings, disorder flags |
| `--mint` | `#43C79F` | Pass states |

**Type:** Archivo (display, used with restraint), Inter (UI), JetBrains Mono
(every PDB ID, sequence, space group, resolution and condition: monospace wherever
data is data).

**Motion:** drawer slide, card flip-out, divider shift, node pulse on cross-highlight.
Nothing else moves. `prefers-reduced-motion` respected throughout.

---

## 2. Data sources (all public, no keys required)

| Source | Endpoint | Used for |
|---|---|---|
| RCSB Search API v2 | `search.rcsb.org/rcsbsearch/v2/query` | Family resolution: sequence (mmseqs2), UniProt accession, Pfam/InterPro annotation, text |
| RCSB Data API | `data.rcsb.org/graphql` | Batch entry / polymer_entity / assembly / chem_comp metadata |
| RCSB 1D Coordinates | `1d-coordinates.rcsb.org/graphql` | UniProt ↔ PDB residue-level mapping |
| RCSB Alignment API | `alignment.rcsb.org/api/v1/` | Pairwise + multiple structure alignment (TM-align, jFATCAT, CE) without local compute |
| RCSB sequence clusters | Data API | Redundancy grouping at 30/50/70/90/95/100 % identity |
| PDBe REST | `ebi.ac.uk/pdbe/api/` | SIFTS, secondary structure, experiment/crystallisation, validation, binding sites |
| PDBe Graph API | `ebi.ac.uk/pdbe/graph-api/` | PDBe-KB annotations, ligand interactions, similar proteins |
| InterPro API | `ebi.ac.uk/interpro/api/` | Pfam, CATH-Gene3D, SUPERFAMILY, PANTHER domain assignments |
| UniProt REST | `rest.uniprot.org/` | Canonical sequence, features, natural variants, isoforms |
| RCSB file service | `files.rcsb.org/download/` | mmCIF for local DSSP / PLIP |
| EMDB | via Data API cross-ref | Map availability, cryo-EM resolution method |
| AlphaFold DB (EBI) | `alphafold.ebi.ac.uk/api/` | Model overlay, pLDDT, disorder cross-check (ties to AlphaFraud) |

**Local compute:** DSSP (secondary structure, topology input), PLIP (interactions),
biotite / gemmi (mmCIF parsing, superposition), parasail or MAFFT (MSA).

**Reuse:** `cif_to_plip.py` already exists and already documents the `pdb_tidy` CONECT
serial-gap bug that corrupts OpenBabel bond perception. Import it, do not rewrite it.

---

## 3. Four-phase rollout

Each phase is independently shippable. Deploy at the end of every phase. The divider
rail carries every phase's section from day one, so nothing needs re-architecting.

---

### PHASE 1 — Spine: resolve a family, map it, ship it

**Goal:** paste anything, get a complete, sortable, exportable inventory of the family,
in the agreed shell. Live on `codswallop.mdeller.com` at the end of this phase.

**Build:**
1. Flask app skeleton: blueprint layout, gunicorn + systemd unit, nginx vhost, certbot.
   Mirror the AlphaFraud deployment pattern exactly.
2. SQLite cache with a `family` / `entry` / `entity` schema and a TTL. Every external
   API response cached; the second query for a family must be instant.
3. **Family resolver** accepting any of: PDB ID, UniProt accession, gene name, Pfam /
   InterPro ID, raw sequence (FASTA or bare), free-text name. Ambiguous input returns a
   disambiguation card, never a guess.
4. Family definition controls: identity threshold slider (30–100 %), include/exclude
   orthologues, include/exclude chimeras and fusion constructs.
5. Batch GraphQL fetch of core metadata for every entry:
   PDB ID · title · method · resolution · R-work · R-free · space group · unit cell ·
   deposition/release date · organism · expression system · chain count · assembly ·
   ligands present · primary citation.
6. **The shell, complete:** drawer front, divider rail with phase pips, family header
   with FILED stamp, stat strip, slide-over entry card, dark/light tokens. Built once,
   here, so later phases only add panels.
7. **Hero map, placeholder layout:** nodes positioned by sequence identity to the family
   centroid until real structural alignment arrives in Phase 3. Entry list beside it,
   cross-highlight wired, Map/Table toggle, method and holo filter chips.
8. Full table view (Tabulator): sortable, filterable, column-toggleable, with
   CSV / JSON / BibTeX export.

**Ships as:** a working family inventory that already beats trawling RCSB by hand.

---

### PHASE 2 — Sequence, construct and domain intelligence

**Goal:** the layer nobody else provides. What was actually *made*, versus what the gene
says.

**Build:**
1. Fetch UniProt canonical; fetch every deposited SEQRES (`entity_poly`).
2. Family MSA (MAFFT or parasail), rendered as an interactive alignment viewer with
   per-column conservation and a sequence logo.
3. **Construct diff engine** — the flagship feature. For every entity, align SEQRES to
   canonical and classify:
   - N-/C-terminal expression tags (His6/His8/His10, Strep-II, FLAG, HA, Myc, Avi)
   - cleavage sites (TEV, 3C/PreScission, thrombin, SUMO, enterokinase) and scar residues
   - fusion partners (MBP, GST, SUMO, Trx, T4 lysozyme, BRIL/apocytochrome b562RIL,
     rubredoxin, PGS, GFP)
   - truncations and internal loop deletions
   - point mutations, sub-classified: catalytic (inactivating), thermostabilising,
     surface entropy reduction (K/E→A patches), disulphide engineering, seleno-Met
   - non-canonical residues and SeMet substitution
4. **Construct table** — one row per unique construct, with the entries that used it and
   their resolutions. Answers "which construct crystallises best?" directly.
5. Domain / fold annotation: Pfam, InterPro, CATH, SCOP, Gene3D per entity via SIFTS +
   InterPro, drawn as a domain-architecture ribbon.
6. **Coverage and disorder census:** stacked per-residue plot across the whole family of
   (a) residues present in construct, (b) residues with modelled density. Highlight
   regions nobody has ever resolved. This feeds the flagged stat in the header.
7. Cross-species / orthologue matrix: which organisms have structures, at what coverage.
8. Wire all of the above into the shared selection state: selecting a construct dims
   every other node on the map.

**Ships as:** the construct-design guide a new project lead actually needs.

---

### PHASE 3 — Structure, chemistry and interactions

**Goal:** see it, superpose it, understand what binds where, and make the hero map real.

**Build:**
1. **Mol\* viewer** embedded, with:
   - single-entry view driven from the slide-over card
   - **family superposition** view (RCSB Alignment API for the transforms, or local
     gemmi/biotite fallback), colour-by-entry or colour-by-conservation
   - toggle: AlphaFold model overlay with pLDDT colouring
2. **Pairwise RMSD / TM-score matrix** across the family, rendered as a clustered
   heatmap with a dendrogram. **This replaces the Phase 1 placeholder layout: node
   positions on the hero map now come from a real structural embedding**, and the
   clusters it identifies (apo vs holo, open vs closed) become first-class filters.
3. **Ligand panel:** every chemical component in the family, with CCD ID, formula,
   SMILES, 2D depiction, occurrence count, and an explicit
   **ligand vs cryoprotectant vs buffer vs ion** classification (PEG, glycerol, sulphate,
   Tris, MES, MPD, cacodylate flagged and separable). Cross-refs to ChEMBL/DrugBank
   where available.
4. **PLIP interactions family-wide:** download mmCIF, convert, run PLIP, store contacts.
   Render as (a) per-entry interaction diagram, (b) a family **interaction fingerprint
   heatmap** (ligand × residue), (c) a "hot residue" ranking of the most frequently
   contacted positions, mapped back onto the alignment and onto Mol\*. Run as a
   background job with progressive fill: the map must not wait for PLIP.
5. **Crystallisation panel:** parse `_exptl_crystal_grow` into structured fields
   (precipitant, salt, buffer, pH, temperature, method, additives). Scatter of pH vs
   precipitant class coloured by resolution; a "what worked" summary table. Export in a
   shape the Top96 predictor can ingest.
6. **Quality / validation panel:** clashscore, RSRZ outliers, Ramachandran, rotamer
   outliers, R-free minus R-work gap, EDS availability, structure-factor and raw-data
   availability. A blunt traffic-light triage of which entries to trust.

**Ships as:** the full scientific payload, and the hero map finally meaning what it says.

---

### PHASE 4 — Topology, export and polish

**Goal:** finish the set, and make it shareable.

**Build:**
1. **Topology diagrams:** generate from DSSP secondary structure as custom SVG
   (pro-origami style: strands as arrows, helices as cylinders, connectivity preserved).
   Build these in-house rather than scraping PDBsum, which is copyrighted.
2. Assembly / oligomeric state panel: author-assigned versus PISA-predicted, with
   disagreements flagged. Interface areas across the family.
3. **Dossier export:** one-click self-contained HTML (and PDF) family report, branded,
   suitable for a project kickoff deck or a grant appendix. Plus BibTeX of every primary
   citation.
4. Shareable permalinks: `codswallop.mdeller.com/f/<family-slug>` with the filter state
   encoded in the URL.
5. Performance: pre-warm popular families, background refresh job, progressive rendering
   so the map and table paint before PLIP finishes.
6. Ship-out: `marcs-vibe-icon` icon (brass-on-slate to match the theme), house-standard
   README, a marcdeller.com blog post, and a link from the mdeller.com landing page.

**Ships as:** the finished, shareable, differentiated tool.

---

## 4. Backlog / stretch (deliberately out of the four phases)

- Conservation mapped onto surface in Mol\* (ConSurf-style, computed in-house)
- Pocket detection (fpocket) and druggability scoring across the family
- Cryo-EM specifics: map resolution vs model resolution, half-map FSC, local resolution
- Nomenclature reconciliation panel (one protein, fourteen names)
- "Who works on this family" — depositing groups, PIs, timeline of labs
- Point-mutation impact overlay from ClinVar / gnomAD for human targets
- Watchlist: email when a new entry lands in a saved family (reuse the AlphaFraud
  weekly timer)
- API endpoint so other tools (BoltzMaker, chatPDB) can consume a family summary as JSON

---

## 5. Ground rules for Claude Code

- British English throughout (colour, licence, behaviour, organise, crystallisation).
- No em dashes: use colons or parentheses.
- Forbidden words: groundbreaking, revolutionary, paradigm-shifting, game-changing.
- Theme tokens exactly as specified in section 1.4. Dark is the default. Never introduce
  a colour outside that table without adding it to the table first.
- Anything physical is brass; all data is cool. That rule is what holds the theme
  together, so do not decorate data panels with brass or dress the drawer in blue.
- The cross-highlight in section 1.3 is load-bearing. Every panel added in later phases
  binds to the same shared selection state rather than keeping its own.
- Mobile-responsive at every phase, not retrofitted. `prefers-reduced-motion` respected.
- Every external link carries `rel="noopener noreferrer"`.
- Cache aggressively and be a polite API citizen: batch GraphQL, respect rate limits,
  set a descriptive User-Agent, back off on 429.
- Attribute data sources visibly (RCSB PDB, PDBe, UniProt, InterPro, PLIP) with
  citations in the footer and the README.
- Never overwrite a working version: version files, keep the previous one.
