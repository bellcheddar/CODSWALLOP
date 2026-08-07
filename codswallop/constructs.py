"""The construct diff engine: what was actually *made*, versus what the gene says.

Aligns every distinct deposited SEQRES against the UniProt canonical and classifies the
differences. This is the layer nobody else provides, and it is the one that answers the
question a project lead actually has: which construct crystallises best, and what did the
people who got it to work put on the ends?

The classification is deliberately conservative. Every motif below is one somebody genuinely
uses in a real vector, and anything unrecognised is reported as an unnamed extension with its
sequence rather than being guessed at or quietly dropped. A wrong name on a fusion partner is
worse than no name: it would send a reader off to clone the wrong thing.
"""

from __future__ import annotations

import re
from typing import Optional

from Bio import Align
from Bio.Align import substitution_matrices

# Bumped whenever the motif tables or the classification change, because the per-family
# construct analysis is cached: without it, a new fusion partner added below would never
# appear for any family already in the cache. Same trap as rcsb.PARSE_VERSION.
ENGINE_VERSION = 4


# --------------------------------------------------------------------------------------
# Reference motifs
# --------------------------------------------------------------------------------------
# Affinity and epitope tags. Ordered longest-first at match time so His10 is not reported as
# His6 with an overhang. The His run is a regex because 6, 8, 10 and the occasional 7 or 9
# all appear, and because vectors interrupt them (HHHHHHSSGHHHHHH in some pET derivatives).
TAGS: list[tuple[str, str]] = [
    ("Strep-II",        "WSHPQFEK"),
    ("FLAG",            "DYKDDDDK"),
    ("HA",              "YPYDVPDYA"),
    ("Myc",             "EQKLISEEDL"),
    ("Avi",             "GLNDIFEAQKIEWHE"),
    ("V5",              "GKPIPNPLLGLDST"),
    ("T7",              "MASMTGGQQMG"),
    ("S-tag",           "KETAAAKFERQHMDS"),
    ("Calmodulin-BP",   "KRRWKKNFIAVSAANRFKKISSSGAL"),
    ("HAT",             "KDHLIHNVHKEFHAHAHNK"),
]
HIS_RUN = re.compile(r"H{5,12}")

# Protease sites, given as the recognition sequence. `scar` is what is left on the protein
# after cleavage, which is the part that actually ends up in the crystal and the part people
# forget about when they wonder why their construct has a stray GP on the front.
PROTEASES: list[tuple[str, str, str]] = [
    ("TEV",           "ENLYFQ",        "G or S"),
    ("3C/PreScission", "LEVLFQ",       "GP"),
    ("Thrombin",      "LVPR",          "GS"),
    ("Factor Xa",     "IEGR",          "none"),
    ("Enterokinase",  "DDDDK",         "none"),
    ("SUMO/Ulp1",     "HSTV",          "none"),
    ("Sortase A",     "LPETG",         "LPET"),
]

# Fusion partners, identified by a distinctive internal peptide rather than the whole
# sequence: the whole sequence would fail on the many point-mutated variants in use
# (BRIL and T4 lysozyme in particular are almost always thermostabilised further).
FUSIONS: list[tuple[str, str, int]] = [
    ("MBP (maltose-binding protein)", "KIEEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDIIFWAHDRFGGYAQSGLLAEITPDKAFQDKLYPFTWDAVRYNGKLIAYPIAVEALSLIYNKDLLPNPPKTWEEIPALDKELKAKGKSALMFNLQEPYFTWPLIAADGGYAFKYENGKYDIKDVGVDNAGAKAGLTFLVDLIKNKHMNADTDYSIAEAAFNKGETAMTINGPWAWSNIDTSKVNYGVTVLPTFKGQPSKPFVGVLSAGINAASPNKELAKEFLENYLLTDEGLEAVNKDKPLGAVALKSYEEELAKDPRIAATMENAQKGEIMPNIPQMSAFWYAVRTAVINAASGRQTVDEALKDAQTNSSSNNNNNNNNNNLGIE", 370),
    ("GST (glutathione S-transferase)", "SPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYIDGDVKLTQSMAIIRYIADKHNMLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLKVDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALDVVLYMDPMCLDAFPKLVCFKKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPK", 210),
    ("SUMO (Smt3)", "SDSEVNQEAKPEVKPEVKPETHINLKVSDGSSEIFFKIKKTTPLRRLMEAFAKRQGKEMDSLRFLYDGIRIQADQTPEDLDMEDNDIIEAHREQIGG", 90),
    ("Thioredoxin (Trx)", "SDKIIHLTDDSFDTDVLKADGAILVDFWAEWCGPCKMIAPILDEIADEYQGKLTVAKLNIDQNPGTAPKYGIRGIPTLLLFKNGEVAATKVGALSKGQLKEFLDANLA", 100),
    ("T4 lysozyme", "NIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL", 155),
    ("BRIL (apocytochrome b562RIL)", "ADLEDNMETLNDNLKVIEKADNAAQVKDALTKMRAAALDAQKATPPKLEDKSPDSPEMKDFRHGFDILVGQIDDALKLANEGKVKEAQAAAEQLKTTRNAYIQKYL", 100),
    ("Rubredoxin", "MQKYVCTVCGYEYDPAEGDPDNGVKPGTSFDDLPADWVCPVCGAPKSEFEAA", 50),
    ("PGS (Pyrococcus glycogen synthase)", "MKIAILGSTGSIGTQTLDVIRHNPDKFKVVGLAAGGNVELLAEQIREFKPKYVAV", 50),
    ("GFP", "SKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK", 230),
    ("Halo tag", "AEIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPMGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISG", 290),
]

# Residues that are not one of the twenty. `MSE` (selenomethionine) is the one that matters
# for phasing and is reported separately; the rest are grouped.
NONCANONICAL = set("BJOUXZ")

# Surface entropy reduction replaces clusters of high-entropy surface residues with alanine.
SER_FROM = set("KEQ")

# Residues whose loss is characteristically a deliberate catalytic knockout when the position
# is annotated as an active site. Used only to phrase the finding, never to invent one.
_STOP = object()


def _aligner() -> Align.PairwiseAligner:
    """Global alignment, BLOSUM62, free terminal gaps.

    Terminal gaps must be free. A construct carrying a 20-residue His-tag and a 40-residue
    truncation is the normal case, and charging for those end gaps makes the aligner shuffle
    the whole alignment to avoid them, which turns one clean truncation into a scatter of
    fictitious internal indels.
    """
    a = Align.PairwiseAligner()
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.mode = "global"
    a.open_gap_score = -11
    a.extend_gap_score = -1
    # The combined setter, deliberately: the separate target_/query_end_gap_score properties
    # are deprecated in Biopython 1.88 and warn on every call.
    a.end_gap_score = 0.0
    return a


_ALIGNER = _aligner()


# --------------------------------------------------------------------------------------
# Motif recognition on the overhangs
# --------------------------------------------------------------------------------------
def _find_tags(segment: str) -> list[dict]:
    """Name every recognisable vector element in one terminal overhang."""
    found: list[dict] = []
    seen_spans: list[tuple[int, int]] = []

    def claim(start: int, end: int) -> bool:
        for s, e in seen_spans:
            if start < e and end > s:
                return False
        seen_spans.append((start, end))
        return True

    for m in HIS_RUN.finditer(segment):
        if claim(m.start(), m.end()):
            found.append({"kind": "tag", "name": f"His{m.end() - m.start()}",
                          "start": m.start(), "end": m.end(), "seq": m.group()})

    for name, motif in TAGS:
        idx = segment.find(motif)
        if idx >= 0 and claim(idx, idx + len(motif)):
            found.append({"kind": "tag", "name": name, "start": idx,
                          "end": idx + len(motif), "seq": motif})

    for name, motif, scar in PROTEASES:
        idx = segment.find(motif)
        if idx >= 0 and claim(idx, idx + len(motif)):
            found.append({"kind": "protease", "name": name, "start": idx,
                          "end": idx + len(motif), "seq": motif, "scar": scar})

    # Fusion partners are matched on a distinctive internal window rather than the full
    # sequence, because the deployed versions are almost always point-mutated.
    for name, seq, min_len in FUSIONS:
        if len(segment) < min_len * 0.55:
            continue
        probe = seq[20:50]
        if probe and probe in segment:
            found.append({"kind": "fusion", "name": name, "start": segment.find(probe),
                          "end": segment.find(probe) + len(probe), "seq": probe})
            continue
        # Fall back to a looser test for the mutated ones: several short windows must hit.
        windows = [seq[i:i + 12] for i in range(0, min(len(seq), 200), 40)]
        hits = sum(1 for w in windows if w and w in segment)
        if hits >= 2:
            found.append({"kind": "fusion", "name": name + " (variant)", "start": 0,
                          "end": len(segment), "seq": None})

    found.sort(key=lambda f: f["start"])
    return found


def _describe_overhang(segment: str, which: str) -> Optional[dict]:
    """Turn a terminal overhang into named parts, keeping whatever is left over visible."""
    if not segment:
        return None
    parts = _find_tags(segment)
    named = sum(p["end"] - p["start"] for p in parts)
    return {
        "terminus": which,
        "length": len(segment),
        "sequence": segment if len(segment) <= 120 else segment[:117] + "...",
        "parts": parts,
        # How much of the overhang we could not account for. A large unexplained overhang is
        # itself worth showing: it is usually an unrecognised fusion or a native extension,
        # and reporting the residues lets a reader recognise it when the engine cannot.
        "unexplained": max(0, len(segment) - named),
    }


# --------------------------------------------------------------------------------------
# The diff
# --------------------------------------------------------------------------------------
def _locate_fusions(construct: str) -> list[dict]:
    """Find fusion partners by searching the construct sequence directly.

    Deliberately independent of the alignment. A 164-residue fusion partner is a 164-residue
    gap to a global aligner, costing about as much as mis-aligning those residues against the
    target, so the aligner routinely shreds a large insertion into several small ones spread
    across the protein: 2RH1's T4 lysozyme came out as a 42-residue insertion after residue
    230 with a scatter of fictitious point mutations either side. Searching for the partner
    itself does not care how the alignment came out.
    """
    hits: list[dict] = []
    taken: list[tuple[int, int]] = []

    def overlaps(s: int, e: int) -> bool:
        return any(s < te and e > ts for ts, te in taken)

    for name, seq, min_len in FUSIONS:
        best: Optional[tuple[int, int]] = None
        # Walk the partner in windows and collect where they land. Windows rather than the
        # whole sequence because the deployed versions are nearly always point-mutated
        # (T4 lysozyme is usually C54T/C97A, BRIL is a thermostabilised b562).
        positions = []
        for i in range(0, len(seq) - 14, 15):
            w = seq[i:i + 15]
            idx = construct.find(w)
            if idx >= 0:
                positions.append((idx, i))
        if len(positions) < 2:
            continue
        # Consistent placement: the windows must appear in the construct in the same order
        # and roughly the same spacing as in the partner, or this is coincidence.
        positions.sort()
        span_start = positions[0][0] - positions[0][1]
        span_end = span_start + len(seq)
        start = max(0, span_start)
        end = min(len(construct), span_end)
        if end - start < min_len * 0.5 or overlaps(start, end):
            continue
        taken.append((start, end))
        hits.append({"kind": "fusion", "name": name, "start": start, "end": end,
                     "matched_windows": len(positions), "length": end - start})
    hits.sort(key=lambda h: h["start"])
    return hits


def diff(canonical: str, construct: str, features: Optional[dict] = None) -> dict:
    """Compare one deposited SEQRES against the canonical sequence.

    `features` is the UniProt feature map (active sites, binding sites, signal peptide), used
    only to *describe* a mutation that lands on one, never to invent a finding.
    """
    features = features or {}
    if not canonical or not construct:
        return {"ok": False}

    # Excise any fusion partner BEFORE aligning, and align only what is left. Leaving it in
    # lets the aligner spread a 164-residue partner across the target as a scatter of fake
    # indels and point mutations; taking it out first gives a clean diff of the protein the
    # family is actually about, plus a named partner and the residue it was spliced after.
    located = _locate_fusions(construct)
    if located:
        kept, cursor, offsets = [], 0, []
        for h in located:
            kept.append(construct[cursor:h["start"]])
            offsets.append((len("".join(kept)), h))
            cursor = h["end"]
        kept.append(construct[cursor:])
        stripped = "".join(kept)
        fusion_sites = offsets
    else:
        stripped, fusion_sites = construct, []

    aln = _ALIGNER.align(canonical, stripped)[0]
    blocks = list(zip(aln.aligned[0], aln.aligned[1]))
    if not blocks:
        return {"ok": False}

    (c_first, _), (q_first, _) = blocks[0]
    (_, c_last), (_, q_last) = blocks[-1]

    result: dict = {
        "ok": True,
        "canonical_length": len(canonical),
        "construct_length": len(construct),
        "stripped_length": len(stripped),
        "canonical_span": [int(c_first) + 1, int(c_last)],
        "n_overhang": _describe_overhang(stripped[:q_first], "N"),
        "c_overhang": _describe_overhang(stripped[q_last:], "C"),
        "mutations": [],
        "deletions": [],
        "insertions": [],
        "noncanonical": [],
        "semet": False,
    }

    # ---- truncations, relative to the canonical ------------------------------------
    if c_first > 0:
        result["deletions"].append({
            "kind": "N-terminal truncation", "start": 1, "end": int(c_first),
            "length": int(c_first),
        })
    if c_last < len(canonical):
        result["deletions"].append({
            "kind": "C-terminal truncation", "start": int(c_last) + 1,
            "end": len(canonical), "length": len(canonical) - int(c_last),
        })

    # ---- internal indels, between consecutive aligned blocks -----------------------
    for i in range(len(blocks) - 1):
        (_, c_end), (_, q_end) = blocks[i]
        (c_start, _), (q_start, _) = blocks[i + 1]
        c_gap, q_gap = int(c_start) - int(c_end), int(q_start) - int(q_end)
        if c_gap > 0:
            result["deletions"].append({
                "kind": "internal deletion", "start": int(c_end) + 1, "end": int(c_start),
                "length": c_gap,
            })
        if q_gap > 0:
            seg = stripped[int(q_end):int(q_start)]
            # Internal insertions get the same motif scan as the terminal overhangs, and
            # this is not a nicety: the canonical GPCR crystallisation trick is to replace
            # intracellular loop 3 with T4 lysozyme or BRIL, so the most famous fusion
            # constructs in the archive carry their partner INTERNALLY. Scanning only the
            # termini reported 2RH1, the original beta-2 adrenergic receptor-T4L structure,
            # as an unfused protein with a handful of point mutations.
            result["insertions"].append({
                "kind": "internal insertion", "after": int(c_end),
                "length": q_gap, "sequence": seg if len(seg) <= 60 else seg[:57] + "...",
                "parts": _find_tags(seg),
            })

    # ---- point mutations, inside aligned blocks ------------------------------------
    act = set(features.get("active_site") or [])
    bind = set(features.get("binding_site") or [])
    for (cs, ce), (qs, qe) in blocks:
        for offset in range(int(ce) - int(cs)):
            c_res = canonical[int(cs) + offset]
            q_res = stripped[int(qs) + offset]
            if c_res == q_res:
                continue
            pos = int(cs) + offset + 1
            mut = {"position": pos, "from": c_res, "to": q_res,
                   "label": f"{c_res}{pos}{q_res}", "classes": []}

            if q_res in NONCANONICAL:
                result["noncanonical"].append(mut["label"])
                continue

            # Catalytic: the position is an annotated active or binding site. Stated as
            # "at an annotated active site", not as "inactivating": whether it killed the
            # enzyme is a claim about an experiment this tool has not seen.
            if pos in act:
                mut["classes"].append("at an annotated active site")
            if pos in bind:
                mut["classes"].append("at an annotated binding site")
            # Surface entropy reduction: a long, flexible, charged residue to alanine.
            if q_res == "A" and c_res in SER_FROM:
                mut["classes"].append("possible surface entropy reduction")
            # Engineered disulphide, or one removed.
            if q_res == "C" and c_res != "C":
                mut["classes"].append("cysteine introduced")
            if c_res == "C" and q_res in "AS":
                mut["classes"].append("cysteine removed")
            # Proline substitutions are the standard prefusion-stabilising trick.
            if q_res == "P" and c_res != "P":
                mut["classes"].append("proline introduced (rigidifying)")
            result["mutations"].append(mut)

    # ---- selenomethionine ----------------------------------------------------------
    # The canonical alphabet has no code for it, so the tell is the entity sequence using
    # the ambiguity codes, or the description saying so. Handled by the caller for the
    # latter; here it is the residue-level test only.
    result["semet"] = "U" in stripped

    result["mutation_count"] = len(result["mutations"])
    result["is_engineered"] = bool(
        result["mutations"] or result["deletions"] or result["insertions"]
        or (result["n_overhang"] and result["n_overhang"]["parts"])
        or (result["c_overhang"] and result["c_overhang"]["parts"])
    )
    # Roll the motif hits up across ALL THREE places a vector element can sit: the two
    # terminal overhangs and any internal insertion.
    def _all_parts():
        for side in ("n_overhang", "c_overhang"):
            if result[side]:
                for p in result[side]["parts"]:
                    yield side[0].upper(), p
        for ins in result["insertions"]:
            for p in ins.get("parts") or []:
                yield "internal", p

    parts = list(_all_parts())

    # The directly located partners, placed relative to the aligned target so the summary can
    # say whether the partner replaced an internal loop or was hung off an end.
    result["fusion_sites"] = []
    for offset, h in fusion_sites:
        if offset <= (q_first or 0) + 2:
            where = "N-terminal"
        elif offset >= q_last - 2:
            where = "C-terminal"
        else:
            where = "internal"
        # Which canonical residue it was spliced after: walk the aligned blocks for the
        # target position matching this offset in the stripped sequence.
        after = None
        for (cs, ce), (qs, qe) in blocks:
            if qs <= offset <= qe:
                after = int(cs) + (offset - int(qs))
                break
        result["fusion_sites"].append({"name": h["name"], "where": where,
                                       "length": h["length"], "after": after})
        parts.append((where, {"kind": "fusion", "name": h["name"]}))
    result["tags"] = sorted({p["name"] for _, p in parts if p["kind"] == "tag"})
    result["proteases"] = sorted({p["name"] for _, p in parts if p["kind"] == "protease"})
    result["fusions"] = sorted({p["name"] for _, p in parts if p["kind"] == "fusion"})
    result["internal_fusions"] = sorted({
        p["name"] for where, p in parts if p["kind"] == "fusion" and where == "internal"})
    result["has_fusion"] = bool(result["fusions"])
    return result


def summarise(d: dict) -> str:
    """One line describing a construct, for the table and the index card."""
    if not d.get("ok"):
        return "not aligned"
    bits = []
    if d["tags"]:
        bits.append("+".join(d["tags"]))
    if d["fusions"]:
        internal_set = set(d.get("internal_fusions") or [])
        names = []
        for f in d["fusions"]:
            short = f.split(" (")[0]
            names.append(short + (" (internal)" if f in internal_set else ""))
        bits.append("fused to " + ", ".join(names))
    if d["proteases"]:
        bits.append(d["proteases"][0] + " site")
    trunc = [x for x in d["deletions"] if "truncation" in x["kind"]]
    if trunc:
        bits.append("residues {}-{}".format(*d["canonical_span"]))
    internal = [x for x in d["deletions"] if x["kind"] == "internal deletion"]
    if internal:
        bits.append(f"{len(internal)} internal deletion" + ("s" if len(internal) > 1 else ""))
    ins = d.get("insertions") or []
    unnamed = [x for x in ins if not (x.get("parts") or [])]
    if unnamed:
        longest = max(unnamed, key=lambda x: x["length"])
        bits.append(f"{longest['length']}-residue insertion after {longest['after']}")
    n = d["mutation_count"]
    if n:
        labels = [m["label"] for m in d["mutations"][:3]]
        bits.append(", ".join(labels) + (f" +{n - 3} more" if n > 3 else ""))
    return "; ".join(bits) if bits else "matches the canonical sequence"
