from __future__ import annotations

"""Coordinate-independent topology grouping for standalone GRO inputs."""

from collections.abc import Iterable
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from string import ascii_uppercase
from typing import Any

from ..models import Atom, Frame
from .trajectory import read_gro


TOPOLOGY_GROUP_LABELS = tuple(ascii_uppercase)
MAX_TOPOLOGY_GROUPS = len(TOPOLOGY_GROUP_LABELS)
TOPOLOGY_FINGERPRINT_SCHEMA = "sqq-gro-topology-v1"


@dataclass(frozen=True)
class GroResidueBlock:
    """One contiguous GRO residue block without its numeric residue id."""

    resname: str
    atomnames: tuple[str, ...]


@dataclass(frozen=True)
class GroTopologyDescriptor:
    """Coordinate-independent atom and molecule layout for one GRO frame."""

    atom_count: int
    residue_blocks: tuple[GroResidueBlock, ...]

    def canonical_payload(self) -> dict[str, Any]:
        """Return the versioned JSON-compatible value used for fingerprinting."""
        return {
            "schema": TOPOLOGY_FINGERPRINT_SCHEMA,
            "atom_count": int(self.atom_count),
            "residue_blocks": [
                {
                    "resname": block.resname,
                    "atomnames": list(block.atomnames),
                }
                for block in self.residue_blocks
            ],
        }


@dataclass(frozen=True)
class GroInputScan:
    """Successful topology pre-scan for one source in finalized input order."""

    source_index: int
    path: Path
    descriptor: GroTopologyDescriptor
    fingerprint: str


@dataclass(frozen=True)
class GroInputFailure:
    """Non-strict pre-scan failure that cannot be assigned to a topology group."""

    source_index: int
    path: Path
    error: str
    exception_type: str


@dataclass(frozen=True)
class GroPrescanResult:
    """Successful and failed records from a read-only GRO pre-scan."""

    inputs: tuple[GroInputScan, ...]
    failures: tuple[GroInputFailure, ...]

    @property
    def source_count(self) -> int:
        return len(self.inputs) + len(self.failures)


@dataclass(frozen=True)
class GroInputAssignment:
    """One successfully scanned source assigned to a first-seen topology group."""

    source_index: int
    path: Path
    descriptor: GroTopologyDescriptor
    fingerprint: str
    group_index: int
    group_label: str | None


@dataclass(frozen=True)
class GroTopologyGroup:
    """All GRO sources sharing one topology fingerprint."""

    group_index: int
    label: str | None
    fingerprint: str
    descriptor: GroTopologyDescriptor
    inputs: tuple[GroInputAssignment, ...]

    @property
    def paths(self) -> tuple[Path, ...]:
        return tuple(item.path for item in self.inputs)

    @property
    def source_indices(self) -> tuple[int, ...]:
        return tuple(item.source_index for item in self.inputs)


@dataclass(frozen=True)
class GroGroupingResult:
    """Stable first-occurrence groups plus source failures and limit metadata."""

    groups: tuple[GroTopologyGroup, ...]
    assignments: tuple[GroInputAssignment, ...]
    failures: tuple[GroInputFailure, ...]
    group_limit: int = MAX_TOPOLOGY_GROUPS
    over_group_limit: bool = False

    @property
    def group_count(self) -> int:
        return len(self.groups)

    @property
    def source_count(self) -> int:
        return len(self.assignments) + len(self.failures)

    @property
    def labels_enabled(self) -> bool:
        """Whether every group has a valid A-Z label."""
        return not self.over_group_limit

    @property
    def info_only_fallback_required(self) -> bool:
        """Whether the caller must apply the whole-run info-only fallback."""
        return self.over_group_limit

    def limit_metadata(self) -> dict[str, int | bool]:
        """Return YAML-safe metadata for root run configuration output."""
        return {
            "topology_group_count": self.group_count,
            "topology_group_limit": self.group_limit,
            "topology_group_limit_exceeded": self.over_group_limit,
            "topology_group_labels_enabled": self.labels_enabled,
            "info_only_fallback_required": self.info_only_fallback_required,
        }

    def source_mapping(self) -> tuple[dict[str, Any], ...]:
        """Return a complete, input-ordered, YAML-safe source mapping."""
        records: list[dict[str, Any]] = [
            {
                "source_index": item.source_index,
                "source": str(item.path),
                "status": "ok",
                "fingerprint": item.fingerprint,
                "group_index": item.group_index,
                "group_label": item.group_label,
            }
            for item in self.assignments
        ]
        records.extend(
            {
                "source_index": item.source_index,
                "source": str(item.path),
                "status": "failed",
                "fingerprint": None,
                "group_index": None,
                "group_label": None,
                "error": item.error,
                "exception_type": item.exception_type,
            }
            for item in self.failures
        )
        records.sort(key=lambda item: int(item["source_index"]))
        return tuple(records)


def topology_group_label(group_index: int) -> str:
    """Return the A-Z label for a zero-based group index."""
    try:
        index = int(group_index)
    except (TypeError, ValueError) as exc:
        raise ValueError("Topology group index must be an integer from 0 through 25.") from exc
    if isinstance(group_index, float) and not group_index.is_integer():
        raise ValueError("Topology group index must be an integer from 0 through 25.")
    if not 0 <= index < MAX_TOPOLOGY_GROUPS:
        raise ValueError(
            f"Topology group index must be between 0 and {MAX_TOPOLOGY_GROUPS - 1}; "
            f"got {group_index!r}."
        )
    return TOPOLOGY_GROUP_LABELS[index]


def gro_topology_descriptor(frame: Frame) -> GroTopologyDescriptor:
    """Describe ordered contiguous residue blocks without coordinate or id values."""
    blocks: list[GroResidueBlock] = []
    current_key: tuple[int, str] | None = None
    current_resname = ""
    current_atomnames: list[str] = []

    for atom in frame.atoms:
        resname = _identity_name(atom.resname, "residue")
        atomname = _identity_name(atom.atomname, "atom")
        # Numeric residue ids identify source block boundaries but are never stored
        # in, or hashed as part of, the topology descriptor.
        key = (int(atom.resid), resname)
        if current_key is not None and key != current_key:
            blocks.append(
                GroResidueBlock(
                    resname=current_resname,
                    atomnames=tuple(current_atomnames),
                )
            )
            current_atomnames = []
        current_key = key
        current_resname = resname
        current_atomnames.append(atomname)

    if current_key is not None:
        blocks.append(
            GroResidueBlock(
                resname=current_resname,
                atomnames=tuple(current_atomnames),
            )
        )

    return GroTopologyDescriptor(
        atom_count=len(frame.atoms),
        residue_blocks=tuple(blocks),
    )


def fingerprint_topology_descriptor(descriptor: GroTopologyDescriptor) -> str:
    """Return a deterministic SHA-256 digest for a GRO topology descriptor."""
    payload = json.dumps(
        descriptor.canonical_payload(),
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("ascii")
    return sha256(payload).hexdigest()


def gro_topology_fingerprint(frame: Frame) -> str:
    """Return the coordinate-independent topology fingerprint for one GRO frame."""
    return fingerprint_topology_descriptor(gro_topology_descriptor(frame))


def scan_gro_inputs(
    paths: Iterable[str | Path],
    *,
    strict: bool = False,
) -> GroPrescanResult:
    """Read GRO sources once in caller order and retain non-strict failures."""
    inputs: list[GroInputScan] = []
    failures: list[GroInputFailure] = []
    for source_index, raw_path in enumerate(paths):
        path = Path(raw_path)
        try:
            if path.suffix.lower() != ".gro":
                raise ValueError(f"GRO topology pre-scan requires a .gro source: {path}")
            frame = read_gro(path)
            descriptor = gro_topology_descriptor(frame)
            fingerprint = fingerprint_topology_descriptor(descriptor)
        except Exception as exc:
            if strict:
                raise
            failures.append(
                GroInputFailure(
                    source_index=source_index,
                    path=path,
                    error=str(exc),
                    exception_type=type(exc).__name__,
                )
            )
            continue
        inputs.append(
            GroInputScan(
                source_index=source_index,
                path=path,
                descriptor=descriptor,
                fingerprint=fingerprint,
            )
        )
    return GroPrescanResult(inputs=tuple(inputs), failures=tuple(failures))


def group_gro_scans(
    scans: Iterable[GroInputScan],
    failures: Iterable[GroInputFailure] = (),
) -> GroGroupingResult:
    """Group successful scans by fingerprint in first-occurrence input order."""
    ordered_scans = tuple(scans)
    ordered_failures = tuple(failures)
    _validate_unique_source_indices(ordered_scans, ordered_failures)

    group_index_by_fingerprint: dict[str, int] = {}
    descriptors: list[GroTopologyDescriptor] = []
    members: list[list[GroInputScan]] = []
    for scan in ordered_scans:
        group_index = group_index_by_fingerprint.get(scan.fingerprint)
        if group_index is None:
            group_index = len(descriptors)
            group_index_by_fingerprint[scan.fingerprint] = group_index
            descriptors.append(scan.descriptor)
            members.append([])
        elif scan.descriptor != descriptors[group_index]:
            raise ValueError(
                "GRO topology fingerprint collision or inconsistent scan record: "
                f"{scan.path}"
            )
        members[group_index].append(scan)

    over_group_limit = len(members) > MAX_TOPOLOGY_GROUPS
    assignments: list[GroInputAssignment] = []
    groups: list[GroTopologyGroup] = []
    for group_index, group_members in enumerate(members):
        # More than 26 groups invokes one whole-run fallback, so no partial A-Z
        # labeling is exposed for only the first 26 groups.
        label = None if over_group_limit else topology_group_label(group_index)
        group_assignments = tuple(
            GroInputAssignment(
                source_index=scan.source_index,
                path=scan.path,
                descriptor=scan.descriptor,
                fingerprint=scan.fingerprint,
                group_index=group_index,
                group_label=label,
            )
            for scan in group_members
        )
        assignments.extend(group_assignments)
        groups.append(
            GroTopologyGroup(
                group_index=group_index,
                label=label,
                fingerprint=group_members[0].fingerprint,
                descriptor=descriptors[group_index],
                inputs=group_assignments,
            )
        )

    assignments.sort(key=lambda item: item.source_index)
    return GroGroupingResult(
        groups=tuple(groups),
        assignments=tuple(assignments),
        failures=tuple(sorted(ordered_failures, key=lambda item: item.source_index)),
        group_limit=MAX_TOPOLOGY_GROUPS,
        over_group_limit=over_group_limit,
    )


def scan_and_group_gro_inputs(
    paths: Iterable[str | Path],
    *,
    strict: bool = False,
) -> GroGroupingResult:
    """Convenience wrapper combining GRO pre-scan and first-seen grouping."""
    prescan = scan_gro_inputs(paths, strict=strict)
    return group_gro_scans(prescan.inputs, prescan.failures)


def _identity_name(value: str, kind: str) -> str:
    text = str(value).strip()
    if not text:
        raise ValueError(f"GRO topology {kind} names must be non-empty.")
    return text.upper()


def _validate_unique_source_indices(
    scans: tuple[GroInputScan, ...],
    failures: tuple[GroInputFailure, ...],
) -> None:
    indexes = [item.source_index for item in scans]
    indexes.extend(item.source_index for item in failures)
    if len(indexes) != len(set(indexes)):
        raise ValueError("GRO topology pre-scan records contain duplicate source indexes.")


__all__ = [
    "MAX_TOPOLOGY_GROUPS",
    "TOPOLOGY_FINGERPRINT_SCHEMA",
    "TOPOLOGY_GROUP_LABELS",
    "GroGroupingResult",
    "GroInputAssignment",
    "GroInputFailure",
    "GroInputScan",
    "GroPrescanResult",
    "GroResidueBlock",
    "GroTopologyDescriptor",
    "GroTopologyGroup",
    "fingerprint_topology_descriptor",
    "gro_topology_descriptor",
    "gro_topology_fingerprint",
    "group_gro_scans",
    "scan_and_group_gro_inputs",
    "scan_gro_inputs",
    "topology_group_label",
]
