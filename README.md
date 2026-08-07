# 🗄️ CODSWALLOP

> **The whole family, filed: every PDB entry for a protein family, cross-referenced in ninety seconds.**

![python](https://img.shields.io/badge/python-3.11-3776AB?logo=python&logoColor=white) ![flask](https://img.shields.io/badge/flask-3.0-000000?logo=flask&logoColor=white) ![phase](https://img.shields.io/badge/phase-1%20of%204-fcb900) ![licence](https://img.shields.io/badge/licence-MIT-467FF7) ![author](https://img.shields.io/badge/author-Marc%20C.%20Deller%2C%20D.Phil.-1C244B)

<table>
<tr>
<td>🌐 <b>Website</b></td><td><a href="https://marcdeller.com" target="_blank" rel="noopener noreferrer">marcdeller.com</a></td>
<td>✉️ <b>Contact</b></td><td><a href="mailto:marc@marcdeller.com">marc@marcdeller.com</a></td>
<td>🐙 <b>GitHub</b></td><td><a href="https://github.com/bellcheddar/CODSWALLOP" target="_blank" rel="noopener noreferrer">bellcheddar/CODSWALLOP</a></td>
</tr>
</table>

---

**C**omparative **O**verview of **D**omains, **S**tructures, **W**hole-families, **A**lignments, **L**igands, **L**ineages, **O**rganisms & **P**DB entries.

A structural biologist has just been handed a protein family and told to get up to speed. CODSWALLOP does the fortnight of manual PDB trawling in ninety seconds: every entry, every construct, every ligand, every crystallisation condition, cross-referenced and interactive.

**Why it matters:** the RCSB and PDBe landing pages are entry-centric. They answer "what is 4XYZ?" very well, and nobody answers the question a working structural biologist actually asks, which is "what do the 87 structures of this family collectively tell me, and what should I do differently?". CODSWALLOP is family-centric: it treats the family as the unit of analysis, so the answer to "which construct crystallises best" or "which residues has nobody ever put in a construct" is one glance rather than a fortnight of tab-juggling. It is useful for: starting a new target, designing a construct, choosing which deposited structure to trust, planning a crystallisation screen, and writing the structural-biology section of a grant.

Live at [codswallop.mdeller.com](https://codswallop.mdeller.com) (once deployed).

---

## 🧭 What makes it different

| # | Differentiator | Why it matters | Phase |
|---|---|---|---|
| 1 | **Construct diffing** | Aligns every deposited SEQRES against the UniProt canonical to expose tags, fusion partners, truncations, surface-entropy-reduction and thermostabilising mutations. This is what actually gets a crystal, and it is buried in mmCIF nobody reads | 2 |
| 2 | **Global cross-filtering** | One selection state shared by every panel. Click a ligand and the map, the entry list, the crystallisation scatter and the coverage plot all filter together | 1 |
| 3 | **Crystallisation intelligence** | Parses `_exptl_crystal_grow` across the whole family into structured conditions, so you can see what actually worked. Feeds the Top96 crystallisation predictor | 3 |
| 4 | **Coverage census** | Which residues has nobody ever seen? A family-wide map answering a construct-design question in one glance | 1 (constructs), 2 (density) |
| 5 | **Interaction fingerprints** | PLIP run family-wide, so binding-site contacts are comparable across entries rather than one-off | 3 |

## ✨ Features (phase 1, live)

- **Ask for anything.** A PDB ID, a PDB entity id, a UniProt accession, a gene name, a Pfam or InterPro accession, a raw sequence (FASTA or bare), or just what the protein is called. Anything genuinely ambiguous returns a disambiguation card rather than a guess: `4HHB` asks whether you meant the alpha or beta globin, `LYZ` asks which of twelve organisms.
- **One definition of a family.** Whatever you type is resolved to a **seed sequence**, and membership is every PDB polymer entity above a percent-identity threshold to that seed. Every member therefore carries a real identity number to the same reference, so the threshold slider means something exact rather than "whatever the annotation happened to say".
- **The drawer.** An archival shell: a divider rail carrying all ten sections from day one, a family header with a rubber-stamped `FILED` mark, a stat strip, and an index card that flips out from any entry anywhere without ever losing the family context.
- **The constellation.** Every entry as a node: outward is decreasing identity to the seed, the sector is the source organism (sized by how many entries it holds), node size is 1/resolution, colour is method, an amber halo means ligand-bound.
- **Cross-highlight.** Hovering any row in the entry list pulses that entry's node brass in the map, and hovering a node highlights its row. One state, two representations, and every panel added in later phases binds to the same shared selection state.
- **Family definition controls.** Percent-identity slider, method chips, and include/exclude toggles for orthologues and for chimeras and fusion constructs. All of it filters client-side and instantly: the family is assembled once at the lowest threshold and every control re-filters the same cached set.
- **Construct coverage census.** How many of the family's constructs contain each residue of the seed, as a depth profile. On lysozyme it finds residues 1 to 18 covered by 66 of 1,686 constructs, which is the signal peptide that is cleaved and never present in the mature protein. On p53 it finds residues 361 to 393 in 10 of 267.
- **Full table view.** Sortable, filterable, column-toggleable, with CSV, JSON and BibTeX export. The BibTeX is deduplicated to one record per paper with the PDB IDs collected in a note, which is the tedious part of assembling a family bibliography by hand.
- **Dark by default, light on request.** A lamp-lit archive theme that honours `prefers-color-scheme` and respects `prefers-reduced-motion`, self-hosted type, and no third-party request on any page load.

## 🧱 Stack

| Layer | Choice | Note |
|---|---|---|
| Backend | Python 3.11, Flask, gunicorn | No compiled toolchain needed on the droplet in phase 1 |
| Cache | SQLite (WAL) | Raw API responses plus assembled families, both with a one-week TTL |
| Front end | Vanilla JS, hand-rolled SVG | The map needs direct control of individual nodes for the cross-highlight pulse |
| Tables | Tabulator 6.3.1 | Vendored, skinned to the theme tokens |
| Type | Archivo, Inter, JetBrains Mono | Self-hosted woff2, latin and latin-ext only, 244 kB total |
| Serving | nginx, certbot | Same droplet pattern as AlphaFraud |

## 📡 Data sources

All public, none require a key. Please cite them, not this tool, when the data does the work.

| Source | Endpoint | Used for |
|---|---|---|
| RCSB Search API v2 | `search.rcsb.org/rcsbsearch/v2/query` | Family resolution: sequence (mmseqs2), UniProt accession, Pfam/InterPro annotation, full text |
| RCSB Data API | `data.rcsb.org/graphql` | Batched entry, polymer entity, assembly and chem_comp metadata |
| UniProt REST | `rest.uniprot.org/` | Canonical sequence, protein name, organism, gene names |
| InterPro API | `ebi.ac.uk/interpro/api/` | Pfam and InterPro accession names and types |

Phases 2 to 4 add the RCSB 1D Coordinates, Alignment and file services, the PDBe REST and Graph APIs, EMDB and the AlphaFold DB.

## 🔧 Installation

```bash
git clone https://github.com/bellcheddar/CODSWALLOP.git
cd CODSWALLOP
python3 -m venv .venv
./.venv/bin/pip install -r requirements.txt
./.venv/bin/python CODSWALLOP.py init
```

## 🚀 Usage

```bash
# Run the development server (defaults to 127.0.0.1:8006)
./.venv/bin/python CODSWALLOP.py serve

# Assemble a family from the shell, without a browser
./.venv/bin/python CODSWALLOP.py build P00918
./.venv/bin/python CODSWALLOP.py build "TEM-1 beta-lactamase"
./.venv/bin/python CODSWALLOP.py build PF00062 --refresh

# Cache housekeeping
./.venv/bin/python CODSWALLOP.py stats
./.venv/bin/python CODSWALLOP.py purge
```

| Command | Does |
|---|---|
| `init` | Create the cache database and apply any column migrations |
| `serve` | Run the development server (`--bind`, `--debug`) |
| `build <query>` | Assemble a family and print a summary. Exits 2 with a candidate list if the query is ambiguous |
| `stats` | Cache statistics: families, entries, entities, cached responses, database size |
| `purge` | Drop expired cache entries |

### Routes

| Route | Does |
|---|---|
| `/` | The closed drawer: one search field |
| `/lookup?q=` | Resolve the input, then redirect to the family or ask which one you meant |
| `/f/<slug>` | The family. Paints complete if cached, paints its filing state and fetches if not |
| `/api/family/<slug>` | Assemble (or serve) a family as JSON |
| `/api/stats` | Cache statistics, and the per-app hit beacon for the mdeller.com launcher |
| `/healthz` | Liveness |

## 🛠️ Configuration

Copy `.env.example` to `.env`. Everything has a sensible default in `codswallop/config.py`.

| Variable | Default | Does |
|---|---|---|
| `DROPLET_SSH` | | SSH destination for `deploy/deploy.sh` |
| `DROPLET_PATH` | `/opt/codswallop` | Where the code lives on the droplet |
| `SERVER_NAME` | `codswallop.mdeller.com` | Hostname certbot issues a certificate for |
| `BIND_ADDR` | `127.0.0.1:8006` | Address gunicorn binds to (8000 to 8005 are taken by the sibling apps) |
| `CACHE_TTL_HOURS` | `168` | How long a cached API response stays fresh |
| `FAMILY_TTL_HOURS` | `168` | How long an assembled family stays fresh |
| `MAX_FAMILY_ENTITIES` | `2500` | Cap on entities pulled into one family |

The two TTLs default to a week because the PDB releases weekly: a family cannot gain a member mid-week.

## 📊 Performance

Measured on carbonic anhydrase II (`P00918`), 1,490 entries, and lysozyme (`P00698`), 1,686 entries.

| Path | Cost |
|---|---|
| Cold family assembly | 7 to 12 s (one sequence search plus 30 batched GraphQL calls) |
| Warm family, served from SQLite | 0.07 s |
| Assembled family payload | ~1.7 MB JSON, ~240 kB gzipped |

The payload is deduplicated before it goes over the wire: sequences, primary citations, chemical components and domain annotations all move to family-level lookup tables, because a family is mostly the same construct solved again. Carbonic anhydrase II has 1,490 entities and 232 distinct sequences. Doing this took the page from 6.3 MB to 1.8 MB.

## 🎨 Theme tokens

Anything physical is brass (drawer front, dividers, plate, pins, stamps); all data is cool. That rule is what holds the theme together.

| Token | Dark | Light | Role |
|---|---|---|---|
| `--bg` | `#131924` | `#EDE7DC` | Room |
| `--bg-2` | `#182031` | `#E4DCCC` | Drawer interior, rail |
| `--sky` | `#0F1522` | `#F4F1E9` | Constellation field |
| `--surface` | `#1F2938` | `#FCFAF5` | Cards, panels |
| `--surface-2` | `#273346` | `#F1ECE0` | Raised, panel headers |
| `--line` | `#2E3A4E` | `#D3C7B0` | Borders |
| `--text` | `#E4EAF6` | `#1D2430` | Type |
| `--dim` | `#93A2BC` | `#55627A` | Secondary type |
| `--mute` | `#65748E` | `#7C8AA3` | Tertiary type |
| `--brass` | `#C9A063` | `#8A6D40` | Anything physical |
| `--brass-dk` | `#8A6D40` | `#6B5330` | Anything physical, shaded |
| `--on-brass` | `#1A1408` | `#FDF8EE` | Ink sitting on brass |
| `--accent` | `#5B8CFF` | `#1F5FD8` | Brand blue |
| `--cyan` | `#4FD1E3` | `#0F8FA6` | X-ray |
| `--amber` | `#F5B93F` | `#B47908` | Cryo-EM, ligand halo |
| `--violet` | `#9B7BE8` | `#6B47C4` | NMR |
| `--oxide` | `#E0685C` | `#B8382C` | FILED stamp, warnings, disorder flags |
| `--mint` | `#43C79F` | `#1B8763` | Neutron, pass states |

## 🌐 Deployment

Mirrors the AlphaFraud pattern exactly.

```bash
bash deploy/deploy.sh                                    # from your Mac, after the first provision
sudo SERVER_NAME=codswallop.mdeller.com bash /opt/codswallop/deploy/provision.sh   # once, on the droplet
```

`provision.sh` is idempotent: system packages, service user, venv, systemd unit, nginx site, Let's Encrypt certificate, and the HTTP/2 patch that certbot does not apply on nginx 1.24.

## ✅ To Do

Roadmap for CODSWALLOP, in dependency order: each phase is independently shippable, and the divider rail carries every phase's section from day one so navigation never needs rebuilding. Suggestions welcome.

### Phase 1: spine (shipped)

- [x] **Flask skeleton, deploy pattern and port allocation.** Blueprint layout, gunicorn config, systemd unit, nginx vhost with the post-certbot HTTP/2 patch and the shared long-cache snippet. Port 8006, because 8000 to 8005 are taken by AlphaFraud, ChemSage, ChatPDB, BoltzMaker, FlexAppeal and PANTS
- [x] **SQLite cache with a family/entry/entity schema and a TTL.** Two layers: raw API responses keyed by the request itself, and assembled families. Second query for a family is 0.07 s against 7 to 12 s cold
- [x] **Cache keys carry a parse version.** Because searches are cached in *parsed* form, adding a field to the parser leaves it absent from every cached row and the consuming code fails silently. Adding the alignment spans without bumping it left the coverage figure reading 100 % and the fusion count reading 0, both of which look like answers rather than missing inputs
- [x] **Column migrations, not just `CREATE TABLE IF NOT EXISTS`.** Which silently does nothing to a table that already exists, so a new column reads back NULL on any deployed instance and the app computes a plausible wrong answer
- [x] **Family resolver accepting seven input kinds.** PDB ID, PDB entity id, UniProt accession, gene name, Pfam/InterPro accession, raw sequence, free text. Ambiguous input returns a disambiguation card, never a guess
- [x] **One definition of family membership.** Every input resolves to a seed sequence and membership is a percent-identity threshold against it, so the slider is exact rather than annotation-dependent. Where the *seed* rather than the input is a choice (which member represents a Pfam family), it is picked by a stated rule and the header says which and why
- [x] **Batched GraphQL fetch of core metadata.** Fifty entities per round trip, each batch cached independently so overlapping families share what they have in common
- [x] **R-factor fallback chain.** `_refine.ls_R_factor_R_work` is the field everyone quotes and a good fraction of the archive leaves it empty, putting the number in `ls_R_factor_obs` or `ls_R_factor_all` instead. Reading only the obvious field showed a blank R-work beside a populated R-free on 82 of 1,473 carbonic anhydrase II X-ray entities
- [x] **The drawer shell, complete.** Divider rail with phase pips, family header with the FILED stamp, stat strip, slide-over index card, dark and light token sets. Built once so later phases only add panels
- [x] **Hero map with a placeholder layout.** Positions come from the server so the map, table and any future export agree on one set of coordinates. The radial axis is identity *rank*, not raw identity: a family whose members are 95 % near-identical piles onto one circle under a linear map, which is exactly what the first version drew
- [x] **Cross-highlight wired.** One shared selection state, and every panel reads it. This is the element the interface is remembered by and the thing later phases must bind to rather than reinventing
- [x] **Construct coverage census as a depth profile, not a binary union.** The union answers "has anyone ever put this residue in a construct", which for any well-studied protein is yes everywhere: 140 of the 1,992 spike entities carry the full 1,273 residues, so the union reports perfect coverage for a protein whose cytoplasmic tail is essentially never studied. The depth profile says the useful thing instead
- [x] **Full table view with CSV, JSON and BibTeX export.** BibTeX deduplicated to one record per paper with the PDB IDs collected in a note
- [x] **Payload deduplication.** Sequences, citations, chemical components and domain annotations moved to family-level lookups: 6.3 MB to 1.8 MB, ~240 kB gzipped
- [x] **Mobile responsive and `prefers-reduced-motion` respected**, at this phase rather than retrofitted
- [x] **Deployed to codswallop.mdeller.com** and listed at the top of `mdeller-landing/apps.json`. The beacon is `/api/stats`, and it needed the footer's "filed so far" line to exist first: the beacon was declared before anything on any page requested it, which would have left the launcher's hit count reading zero forever with nothing to show why. A static asset would have been the easier beacon and a worse one, since `/static/` is served `immutable` and a returning reader never re-requests it

### Phase 2: sequence, construct and domain intelligence

- [ ] **Fetch the UniProt canonical and every deposited SEQRES.** The input to everything else in this phase
- [ ] **Family MSA with per-column conservation and a sequence logo.** MAFFT or parasail, rendered as an interactive alignment viewer
- [ ] **Construct diff engine.** The flagship. Align every SEQRES to the canonical and classify: N- and C-terminal expression tags (His6/8/10, Strep-II, FLAG, HA, Myc, Avi), cleavage sites and scar residues (TEV, 3C/PreScission, thrombin, SUMO, enterokinase), fusion partners (MBP, GST, SUMO, Trx, T4 lysozyme, BRIL, rubredoxin, PGS, GFP), truncations and internal loop deletions, point mutations sub-classified as catalytic, thermostabilising, surface entropy reduction or disulphide engineering, and non-canonical residues including SeMet
- [ ] **Replace the phase 1 fusion heuristic.** Currently a blunt length test: a construct 80 or more residues longer than its aligned region is called a fusion. It catches the real ones (7BM4_1 is 332 residues with 252 aligned) and cannot name the partner
- [ ] **Construct table.** One row per unique construct with the entries that used it and their resolutions, answering "which construct crystallises best?" directly
- [ ] **Domain and fold annotation per entity.** Pfam, InterPro, CATH, SCOP and Gene3D via SIFTS and InterPro, drawn as an architecture ribbon
- [ ] **Density coverage census.** The harder question layered onto the same axis as phase 1's construct census: of the residues that *were* in a construct, which has nobody ever seen density for. The two are different questions and the UI must keep labelling them as such
- [ ] **Cross-species orthologue matrix.** Which organisms have structures, at what coverage
- [ ] **Bind it all to the shared selection state.** Selecting a construct dims every other node on the map

### Phase 3: structure, chemistry and interactions

- [ ] **Mol\* viewer embedded.** Single-entry from the index card, plus a family superposition coloured by entry or by conservation, and an AlphaFold overlay with pLDDT colouring. FlexAppeal already vendors Mol\*: reuse it
- [ ] **Pairwise TM-score matrix as a clustered heatmap with a dendrogram.** This retires the placeholder layout: node positions move to a real structural embedding and the clusters it finds (apo versus holo, open versus closed) become first-class filters
- [ ] **Ligand panel.** CCD ID, formula, SMILES, 2D depiction, occurrence count, and an explicit ligand versus cryoprotectant versus buffer versus ion classification. Phase 1 ships a blunt 56-entry exclusion list so the amber halo is roughly right on day one; this replaces it with a per-component verdict
- [ ] **PLIP family-wide.** Per-entry interaction diagrams, a family interaction fingerprint heatmap (ligand by residue), and a hot-residue ranking mapped back onto the alignment and onto Mol\*. Runs as a background job with progressive fill: the map must not wait for PLIP. Reuse `cif_to_plip.py`, which already documents the `pdb_tidy` CONECT serial-gap bug that corrupts OpenBabel bond perception, and do not rewrite it
- [ ] **Crystallisation panel.** Parse `_exptl_crystal_grow` into structured precipitant, salt, buffer, pH, temperature, method and additives. Scatter of pH against precipitant class coloured by resolution, and a "what worked" table exported in a shape TopPDBLX can ingest
- [ ] **Quality and validation panel.** Clashscore, RSRZ outliers, Ramachandran and rotamer outliers, R-free minus R-work gap, EDS and structure-factor availability. A blunt traffic-light triage of which entries to trust

### Phase 4: topology, export and polish

- [ ] **Topology diagrams from DSSP as custom SVG.** Strands as arrows, helices as cylinders, connectivity preserved. Built in-house rather than scraped from PDBsum, which is copyrighted
- [ ] **Assembly and oligomeric state panel.** Author-assigned versus PISA-predicted, with disagreements flagged, and interface areas across the family
- [ ] **Dossier export.** One-click self-contained HTML and PDF family report, branded, suitable for a project kickoff or a grant appendix, plus BibTeX of every primary citation
- [ ] **Shareable permalinks with the filter state encoded in the URL.** `/f/<slug>` already exists; this adds the state
- [ ] **Performance.** Pre-warm popular families, background refresh job, progressive rendering
- [ ] **Ship-out.** Icon, blog post on marcdeller.com, and a link from the mdeller.com landing page

### Backlog (deliberately outside the four phases)

- [ ] **Conservation mapped onto the surface in Mol\***, ConSurf-style, computed in-house
- [ ] **Pocket detection (fpocket) and druggability scoring** across the family
- [ ] **Cryo-EM specifics.** Map resolution versus model resolution, half-map FSC, local resolution
- [ ] **Nomenclature reconciliation panel.** One protein, fourteen names
- [ ] **"Who works on this family".** Depositing groups, PIs, a timeline of labs
- [ ] **Point-mutation impact overlay** from ClinVar and gnomAD for human targets
- [ ] **Watchlist.** Email when a new entry lands in a saved family, reusing the AlphaFraud weekly timer
- [ ] **JSON API** so BoltzMaker and chatPDB can consume a family summary

## 📚 Citing the data

CODSWALLOP is a lens over other people's archives. If the data does the work in your paper, cite them:

- **RCSB PDB**: Burley, S.K. et al. *Nucleic Acids Res.* (2023) `10.1093/nar/gkac1077`
- **PDBe**: Armstrong, D.R. et al. *Nucleic Acids Res.* (2020) `10.1093/nar/gkz990`
- **UniProt**: The UniProt Consortium. *Nucleic Acids Res.* (2023) `10.1093/nar/gkac1052`
- **InterPro**: Paysan-Lafosse, T. et al. *Nucleic Acids Res.* (2023) `10.1093/nar/gkac993`

## 📄 Licence

MIT. See [LICENSE](LICENSE).

---

## 👤 Author

**Marc C. Deller, D.Phil.**  
Structural biologist & drug discovery scientist  

<table>
<tr>
<td>🌐</td><td><a href="https://marcdeller.com" target="_blank" rel="noopener noreferrer">marcdeller.com</a></td>
<td>✉️</td><td><a href="mailto:marc@marcdeller.com">marc@marcdeller.com</a></td>
<td>🐙</td><td><a href="https://github.com/bellcheddar/CODSWALLOP" target="_blank" rel="noopener noreferrer">github.com/bellcheddar/CODSWALLOP</a></td>
</tr>
</table>
