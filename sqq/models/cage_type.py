"""Canonical cage-type names and face-composition parsing."""

from __future__ import annotations


KNOWN_CAGE_TYPES = ["512", "51262", "51263", "51264", "51268", "435663"]

TARGET_FACE_COUNTS = {
    "512": {4: 0, 5: 12, 6: 0},
    "51262": {4: 0, 5: 12, 6: 2},
    "51263": {4: 0, 5: 12, 6: 3},
    "51264": {4: 0, 5: 12, 6: 4},
    "51268": {4: 0, 5: 12, 6: 8},
    "435663": {4: 3, 5: 6, 6: 3},
}

CAGE_REPORT_GROUPS = {
    "I": ("512", "51262"),
    "II": ("512", "51264"),
    "H": ("512", "51268", "435663"),
    "HS-I": ("512", "51262", "51263"),
    "TS-I": ("512", "51262", "51263"),
    "I2II": ("51263",),
}


def parse_cage_face_label(label: str) -> dict[int, int] | None:
    """Parse a named cage, numeric 1-10-2, or generic 4^1-5^10-6^2 label."""
    text = label.strip()
    if not text:
        return None
    if text in TARGET_FACE_COUNTS:
        return dict(TARGET_FACE_COUNTS[text])
    if text.count("-") == 2 and all(part.strip().isdigit() for part in text.split("-")):
        n4, n5, n6 = (int(part) for part in text.split("-"))
        return {4: n4, 5: n5, 6: n6}
    counts: dict[int, int] = {}
    for token in text.replace("_", "-").split("-"):
        token = token.strip()
        if not token:
            continue
        if "^" not in token:
            return None
        size_text, count_text = token.split("^", 1)
        if not size_text.isdigit() or not count_text.isdigit():
            return None
        counts[int(size_text)] = int(count_text)
    return counts or None


def cage_type_for_counts(counts: dict[int, int]) -> str:
    """Use a compact named label when available, otherwise a generic label."""
    for name in KNOWN_CAGE_TYPES:
        if _counts_match(counts, TARGET_FACE_COUNTS[name]):
            return name
    return canonical_cage_face_label(counts)


def canonical_cage_type(label: str) -> str:
    """Normalize one report label to the cage type emitted by the search."""
    text = str(label).strip()
    counts = parse_cage_face_label(text)
    if counts is None:
        raise ValueError(
            f"Unsupported cage type '{label}'. Use a named type such as 51268 "
            "or a face-count label such as 4^1-5^10-6^2."
        )
    return cage_type_for_counts(counts)


def canonical_cage_face_label(counts: dict[int, int]) -> str:
    """Return an ASCII cage type label safe for summary columns and filenames."""
    return "-".join(
        f"{size}^{counts.get(size, 0)}"
        for size in sorted(counts)
        if counts.get(size, 0) > 0
    )


def _counts_match(counts: dict[int, int], target_counts: dict[int, int]) -> bool:
    """Compare face counts while treating missing sizes as zero."""
    sizes = set(counts) | set(target_counts)
    return all(counts.get(size, 0) == target_counts.get(size, 0) for size in sizes)


__all__ = [
    "CAGE_REPORT_GROUPS",
    "KNOWN_CAGE_TYPES",
    "TARGET_FACE_COUNTS",
    "cage_type_for_counts",
    "canonical_cage_face_label",
    "canonical_cage_type",
    "parse_cage_face_label",
]
