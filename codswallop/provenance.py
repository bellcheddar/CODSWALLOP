"""Where a family came from: the names it is deposited under, and the people who did it.

Two questions no entry page can answer, because both are properties of the whole body of
work rather than of any one structure.

**Names.** A protein is deposited under whatever the depositor typed, and over thirty years
that drifts: hen lysozyme arrives as `LYSOZYME`, `Lysozyme C`, `LYSOZYME C`, `Hen egg white
lysozyme` and a dozen more. Searching any one of those finds a fraction of the family. This
panel reconciles every spelling actually used against UniProt's own recommended name, its
alternative names and its gene names, and says which are recognised and which are nobody's
official name at all.

**People.** Structural biology is done by groups, and a family's archive is a record of
which ones. The last author is used as the group's proxy, which is a convention of the
field and not a fact in the data: it is labelled as an inference everywhere it is shown,
and the first author is carried beside it so a reader can see the person who did the work
as well as the person whose name the group goes by.

Both are computed from data already fetched for other panels: entity descriptions, the
citation records the BibTeX export already deduplicates, and deposit dates. No new API.
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Optional

PARSE_VERSION = 1

# Names below this share of the family are pooled rather than listed. A thirty-year archive
# has a long tail of one-off spellings and the panel is about the shape of the drift, not a
# census of typos.
MIN_NAME_SHARE = 0.005
MAX_NAMES = 24
# Groups shown before the tail is pooled.
MAX_GROUPS = 20


def _norm(text: Optional[str]) -> str:
    """Fold a deposited name to something comparable.

    Case, punctuation and whitespace only. Deliberately not stemming or dropping words: the
    difference between "Lysozyme C" and "Lysozyme" is the kind of thing this panel exists to
    show, and normalising it away would be answering the question by deleting it.
    """
    t = (text or "").strip().lower()
    t = re.sub(r"[(),.;:]", " ", t)
    t = re.sub(r"[\s_/-]+", " ", t)
    return t.strip()


def build_names(members: list[dict], seed: dict) -> dict:
    """The names this family is deposited under, in both directions.

    `seed` is the resolved seed record, carrying UniProt's recommended name, any alternative
    names and the gene names.

    **Partitioned by accession first, and this is the whole correctness of the panel.** A
    family is assembled at 30 % identity, so it legitimately contains other proteins:
    counting description strings across all of it reported `ALPHA-LACTALBUMIN` and
    `HUMAN LYSOZYME` as unrecognised names for hen lysozyme, when they are the correct names
    of different proteins that happen to be relatives. Only entities carrying the seed's own
    accession can say anything about what the seed is called.

    Two questions, and the second is the one that actually breaks a literature search:

    * **One protein, many names.** P00698 arrives as `Lysozyme C`, `Lysozyme`, `LYSOZYME C`,
      `LYSOZYME`, `HEN EGG WHITE LYSOZYME`, `PROTEIN (LYSOZYME)` and more.
    * **One name, many proteins.** `Lysozyme C` in this same family is *also* the deposited
      description of P61626, which is human lysozyme. Searching the string finds two
      proteins and no way to tell which is which.
    """
    official: dict[str, str] = {}          # normalised -> what kind of official name it is
    display: dict[str, str] = {}           # normalised -> the official spelling

    def register(value, kind):
        key = _norm(value)
        if not key:
            return
        # First registration wins, and the kinds are registered in order of authority, so a
        # string that is both the recommended name and a gene synonym reads as the former.
        official.setdefault(key, kind)
        display.setdefault(key, value)

    register(seed.get("name"), "recommended")
    for alt in seed.get("alt_names") or []:
        register(alt, "alternative")
    for short in seed.get("short_names") or []:
        register(short, "short")
    for gene in seed.get("genes") or []:
        register(gene, "gene")

    accession = (seed.get("accession") or "").upper()
    def classify(key: str, raw: str) -> Optional[str]:
        """What kind of official name this spelling is, if any.

        Exact match first, then the two qualifier forms UniProt itself uses, because a
        depositor writing `Isoform 2B of GTPase KRas` or `GTPase KRas, N-terminally
        processed` has used the official name and added a true qualifier: they are not the
        drift this panel is looking for. Left unhandled they were a third of KRAS's own
        entities, which would have read as a third of the family calling it something made
        up.

        The qualifiers are stripped from the RAW description and not from the normalised
        key, because `_norm` has already turned every comma into a space by then: the
        comma-qualifier branch written against the key could never once have matched.
        """
        direct = official.get(key)
        if direct:
            return direct
        isoform = re.sub(r"^isoform\s+.+?\s+of\s+", "", (raw or "").strip(),
                         flags=re.IGNORECASE)
        if isoform and _norm(isoform) != key and official.get(_norm(isoform)):
            return "variant"
        head = (raw or "").split(",", 1)[0].strip()
        if head and _norm(head) != key and official.get(_norm(head)):
            return "variant"
        return None

    counts: Counter = Counter()
    spellings: dict[str, Counter] = {}
    # name -> accession -> count, for the collision half of the panel.
    by_name_accession: dict[str, Counter] = {}
    other_accessions: dict[str, Counter] = {}
    n_other = 0

    for m in members:
        raw = (m.get("description") or "").strip()
        acc = (m.get("uniprot") or "").upper()
        if raw:
            key = _norm(raw)
            if acc:
                by_name_accession.setdefault(key, Counter())[acc] += 1
        if accession and acc and acc != accession:
            n_other += 1
            if raw:
                other_accessions.setdefault(acc, Counter())[raw] += 1
            continue
        if not raw:
            continue
        # Entities with no accession at all are counted with the seed's: they are usually
        # the same protein deposited before the cross-reference existed, and dropping them
        # would understate the drift this panel is measuring.
        key = _norm(raw)
        counts[key] += 1
        spellings.setdefault(key, Counter())[raw] += 1

    # One name, many proteins. Only collisions inside this family are reported, because that
    # is the set a reader is actually looking at.
    collisions = []
    for key, accs in by_name_accession.items():
        if len(accs) < 2:
            continue
        best = spellings.get(key)
        label = (best.most_common(1)[0][0] if best
                 else max(accs, key=lambda a: accs[a]))
        if not best:
            # A name used only by other accessions still collides among them.
            label = key
        collisions.append({
            "name": label,
            "total": sum(accs.values()),
            "accessions": [{"accession": a, "n": n} for a, n in accs.most_common()],
        })
    collisions.sort(key=lambda c: -c["total"])

    total = sum(counts.values())
    if not total:
        return {"n": 0, "total": 0, "rows": [], "official": [], "recognised": 0,
                "unrecognised": 0, "pooled": 0, "pooled_names": 0,
                "accession": accession or None, "collisions": collisions[:12],
                "n_collisions": len(collisions), "n_other": n_other, "others": []}

    rows = []
    pooled = pooled_names = 0
    for key, n in counts.most_common():
        # The commonest spelling of this name is what gets shown; the rest are counted.
        variants = spellings[key]
        best = variants.most_common(1)[0][0]
        if n / total < MIN_NAME_SHARE and len(rows) >= 6:
            pooled += n
            pooled_names += 1
            continue
        if len(rows) >= MAX_NAMES:
            pooled += n
            pooled_names += 1
            continue
        rows.append({
            "name": best,
            "n": n,
            "share": round(n / total, 4),
            "kind": classify(key, best),
            "variants": len(variants),
        })

    # Counted over every spelling, not only the listed ones, so the headline figure does not
    # move when the tail is pooled for display.
    recognised = sum(n for key, n in counts.items()
                     if classify(key, spellings[key].most_common(1)[0][0]))
    return {
        "n": len(counts),
        "total": total,
        "rows": rows,
        "pooled": pooled,
        "pooled_names": pooled_names,
        "recognised": recognised,
        "unrecognised": total - recognised,
        # What UniProt says it should be called, so the panel can show the target as well as
        # the drift. Deduplicated on the normalised form, keeping the authoritative spelling.
        "official": [{"name": display[k], "kind": v} for k, v in official.items()],
        "accession": accession or None,
        # One name, many proteins.
        "collisions": collisions[:12],
        "n_collisions": len(collisions),
        # The relatives a 30 % identity search legitimately brought in. Reported as other
        # proteins rather than as misnamed copies of the seed.
        "n_other": n_other,
        "others": [
            {"accession": acc, "n": sum(c.values()), "name": c.most_common(1)[0][0]}
            for acc, c in sorted(other_accessions.items(),
                                 key=lambda kv: -sum(kv[1].values()))[:8]
        ],
    }


def _year(date: Optional[str]) -> Optional[int]:
    if not date:
        return None
    m = re.match(r"(\d{4})", str(date))
    return int(m.group(1)) if m else None


def build_people(entries: list[dict]) -> dict:
    """Who deposited this family, and when.

    Counted per ENTRY rather than per paper. One paper routinely covers a series of
    depositions, so counting papers would say a group that solved thirty structures and
    published them together did less work than one that published three papers about one
    structure each. Both numbers are reported, because the ratio is itself informative.
    """
    groups: dict[str, dict] = {}
    by_year: Counter = Counter()
    n_cited = 0

    for e in entries:
        year = _year(e.get("deposit_date")) or _year(e.get("release_date"))
        if year:
            by_year[year] += 1
        cit = e.get("citation") or {}
        authors = [a for a in (cit.get("authors") or []) if a]
        if not authors:
            continue
        n_cited += 1
        # Last author as the group's proxy. A convention of the field, not a fact in the
        # data, and labelled as such wherever this is displayed.
        pi = authors[-1]
        g = groups.setdefault(pi, {
            "pi": pi, "entries": 0, "papers": set(), "years": [],
            "first_authors": Counter(), "methods": Counter(), "best": None,
        })
        g["entries"] += 1
        key = cit.get("doi") or cit.get("title")
        if key:
            g["papers"].add(key)
        if year:
            g["years"].append(year)
        if authors[0] != pi:
            g["first_authors"][authors[0]] += 1
        if e.get("method"):
            g["methods"][e["method"]] += 1
        res = e.get("resolution")
        if res and (g["best"] is None or res < g["best"]):
            g["best"] = res

    rows = []
    for g in groups.values():
        years = sorted(g["years"])
        rows.append({
            "pi": g["pi"],
            "entries": g["entries"],
            "papers": len(g["papers"]),
            "first": years[0] if years else None,
            "last": years[-1] if years else None,
            "best_resolution": g["best"],
            "top_first_author": (g["first_authors"].most_common(1)[0][0]
                                 if g["first_authors"] else None),
            "method": (g["methods"].most_common(1)[0][0] if g["methods"] else None),
        })
    rows.sort(key=lambda r: (-r["entries"], r["pi"]))

    shown = rows[:MAX_GROUPS]
    tail = rows[MAX_GROUPS:]
    total_entries = len(entries)
    return {
        "n_groups": len(rows),
        "n_entries": total_entries,
        "n_cited": n_cited,
        # Entries with no primary citation at all: they exist, and a panel that silently
        # dropped them would overstate how completely the groups below account for the
        # family.
        "n_uncited": total_entries - n_cited,
        "rows": shown,
        "tail_groups": len(tail),
        "tail_entries": sum(r["entries"] for r in tail),
        # Concentration: how much of the family the top group and the top five account for.
        # This is the number that makes the panel worth having, and it is a family-only
        # question by construction.
        "top_share": round(rows[0]["entries"] / total_entries, 4) if rows and total_entries else 0,
        "top5_share": (round(sum(r["entries"] for r in rows[:5]) / total_entries, 4)
                       if rows and total_entries else 0),
        "timeline": [{"year": y, "n": by_year[y]} for y in sorted(by_year)],
        "span": [min(by_year), max(by_year)] if by_year else None,
    }


def build(fam: dict, members: list[dict], seed: dict) -> dict:
    return {
        "version": PARSE_VERSION,
        "names": build_names(members, seed),
        "people": build_people(fam.get("entries") or []),
    }
