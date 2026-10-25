"""Simultaneous semantic main-figure mapping; no manuscript text is bundled.

Only spans introduced by Figure/Fig. are candidates. Supplementary S numbers
are deliberately outside the grammar. Every mapping returns a reason for an
author-review ledger instead of relying on sequential global replacement.
"""
from dataclasses import dataclass
import re

PREFIX = re.compile(r"\b(?:Figures?|Figs?\.)\s*", re.I)
ITEM = re.compile(r"(?:(?P<number>[1-6])(?:\s*(?P<panel>[A-G]))?|(?P<implicit>[A-G]))(?![A-Za-z0-9])")
SEPARATOR = re.compile(r"\s*(?P<sep>[-–—,]|and|&)\s*")


@dataclass(frozen=True)
class Target:
    figure: int
    panel: str = ""

    def label(self):
        return f"{self.figure}{self.panel}"


def target(old):
    n, panel = old.figure, old.panel
    if n == 1:
        if panel == "A":
            return Target(1), "deleted broad workflow replaced by new overview"
        if panel:
            return Target(2, chr(ord(panel) - 1)), "old 1B–E data panel shifted to 2A–D"
        return Target(2), "old Figure 1 data figure shifted to Figure 2"
    if n == 5:
        if panel == "A":
            return Target(1, "E"), "deleted experimental scheme represented in new module 1E"
        if panel:
            return Target(6, chr(ord(panel) - 1)), "old 5B–G data panel shifted to 6A–F"
        return Target(6), "old Figure 5 experimental data shifted to Figure 6"
    if n in (2, 3, 4):
        return Target(n + 1, panel), "main figure number shifted by one"
    raise ValueError(f"Unexpected old main figure {old.label()}")


def _parse_after(text, pos):
    first = ITEM.match(text, pos)
    if not first or not first.group("number"):
        return None
    current = Target(int(first.group("number")), first.group("panel") or "")
    items, cursor = [current], first.end()
    while True:
        sep = SEPARATOR.match(text, cursor)
        if not sep:
            break
        nxt = ITEM.match(text, sep.end())
        if not nxt:
            break
        if nxt.group("number"):
            item = Target(int(nxt.group("number")), nxt.group("panel") or "")
        else:
            if not current.panel:
                break
            item = Target(current.figure, nxt.group("implicit"))
        separator = sep.group("sep")
        if separator in ("-", "–", "—"):
            if current.figure == item.figure and current.panel and item.panel:
                if ord(item.panel) < ord(current.panel):
                    raise ValueError("Descending panel range")
                items.extend(Target(current.figure, chr(a)) for a in range(ord(current.panel)+1, ord(item.panel)+1))
            elif not current.panel and not item.panel:
                if item.figure < current.figure:
                    raise ValueError("Descending figure range")
                items.extend(Target(a) for a in range(current.figure+1, item.figure+1))
            else:
                raise ValueError("Ambiguous cross-figure range; explicit review required")
        else:
            items.append(item)
        current, cursor = item, nxt.end()
    return cursor, items


def extract(text):
    refs = []
    for prefix in PREFIX.finditer(text):
        parsed = _parse_after(text, prefix.end())
        if parsed:
            end, items = parsed
            refs.append({"start": prefix.start(), "end": end,
                         "text": text[prefix.start():end], "items": items})
    return refs


def _compact(targets):
    """Compact adjacent panels within each figure, preserving encounter order."""
    unique = list(dict.fromkeys(targets))
    # A whole-figure overview encompasses any specific module in the same span.
    whole = {t.figure for t in unique if not t.panel}
    unique = [t for t in unique if not t.panel or t.figure not in whole]
    groups = []
    for item in unique:
        if groups and groups[-1][-1].figure == item.figure and item.panel and groups[-1][-1].panel:
            groups[-1].append(item)
        else:
            groups.append([item])
    rendered = []
    for group in groups:
        n = group[0].figure
        ps = [x.panel for x in group]
        consecutive = len(ps) > 1 and all(ord(b) == ord(a)+1 for a,b in zip(ps,ps[1:]))
        if consecutive:
            rendered.append(f"{n}{ps[0]}–{ps[-1]}")
        elif len(ps) > 1:
            rendered.append(str(n)+", ".join(ps[:-1])+" and "+ps[-1])
        else:
            rendered.append(group[0].label())
    prefix = "Figures " if len(groups) > 1 else "Figure "
    return prefix + (" and ".join(rendered) if len(rendered)<3 else ", ".join(rendered[:-1])+", and "+rendered[-1])


def rewrite(text):
    refs, out, end, rows = extract(text), [], 0, []
    for ref in refs:
        replacements, reasons = zip(*(target(item) for item in ref["items"]))
        new = _compact(replacements)
        out.extend((text[end:ref["start"]], new))
        rows.append({"old_ref": ref["text"], "new_ref": new,
                     "reason": "; ".join(dict.fromkeys(reasons)),
                     "old_targets": ",".join(x.label() for x in ref["items"]),
                     "new_targets": ",".join(x.label() for x in replacements)})
        end = ref["end"]
    out.append(text[end:])
    return "".join(out), rows


def scientific_numbers(text):
    """A conservative numeric guard excluding explicit figure and citation labels."""
    clean = text
    spans = extract(text)
    for item in reversed(spans):
        clean = clean[:item["start"]] + "FIGREF" + clean[item["end"]:]
    clean = re.sub(r"\b(?:Supplementary\s+)?(?:Figures?|Figs?\.)\s*S\d+[A-Z]?(?:\s*[-–,]\s*S?\d*[A-Z]?)?", "SUPPREF", clean, flags=re.I)
    clean = re.sub(r"\[(?:\d+[\s,;–-]*)+\]", "CITATION", clean)
    return re.findall(r"(?<![A-Za-z])[-+−]?\d+(?:[.,]\d+)*(?:[eE][-+]?\d+)?%?", clean)
