#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


DEFAULT_REGIONS = Path("data/inputs/luciferase_hesc_h1_regions.tsv")
DEFAULT_H1_ASSOCIATIONS = Path("outputs/cache/hesc_h1/ld_r2_equal_higher_0.8_hESC_H1_associated_regions.tsv")
DEFAULT_LD_SOURCE = Path("data/inputs/ld_r2_equal_higher_0.8")
DEFAULT_TAG_SOURCE = Path("data/inputs/TECAC_GWAS_index_SNPs_OCT_2025")
DEFAULT_RESOLVED_SNPS = Path("data/inputs/luciferase_hesc_h1_resolved_snps.bed")
DEFAULT_PROVENANCE = Path("data/inputs/luciferase_hesc_h1_provenance.json")
DEFAULT_OUTPUT_PREFIX = "outputs/luciferase_hesc_h1/luciferase_hesc_h1"
DEFAULT_DIRECT_BEDPE = (
    Path("data/public/encode_hesc_h1_loops/ENCFF324UIT.bedpe.gz"),
    Path("data/public/encode_hesc_h1_loops/ENCFF401IWZ.bedpe.gz"),
    Path("data/public/encode_hesc_h1_loops/ENCFF519OAV.bedpe.gz"),
    Path("data/public/encode_hesc_h1_loops/ENCFF753NSM.bedpe.gz"),
)


@dataclass(frozen=True)
class Region:
    region_id: str
    chrom: str
    start_1based: int
    end_1based: int
    snps: tuple[str, ...]
    evidence: str
    note: str

    @property
    def start0(self) -> int:
        return self.start_1based - 1

    @property
    def end0(self) -> int:
        return self.end_1based

    @property
    def region_text(self) -> str:
        return f"{self.chrom}:{self.start_1based}-{self.end_1based}"


@dataclass(frozen=True)
class SnpCoordinate:
    chrom: str
    start0: int
    end0: int


@dataclass(frozen=True)
class H1Link:
    region: Region
    query_variant: str
    target_chrom: str
    target_start0: int
    target_end0: int
    target_region: str
    target_name: str
    method: str
    source_accession: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create focused UCSC Genome Browser tracks for luciferase-active "
            "regions evaluated against verified ENCODE H1 ChIA-PET loops."
        )
    )
    parser.add_argument(
        "--regions",
        type=Path,
        default=DEFAULT_REGIONS,
        help=f"Region TSV with region_id/chrom/start_1based/end_1based/snps columns (default: {DEFAULT_REGIONS})",
    )
    parser.add_argument(
        "--h1-associations",
        type=Path,
        default=DEFAULT_H1_ASSOCIATIONS,
        help=(
            "Cached H1-hESC associated-region TSV to filter by requested SNPs "
            f"(default: {DEFAULT_H1_ASSOCIATIONS})"
        ),
    )
    parser.add_argument(
        "--source-mode",
        choices=["direct", "cached"],
        default="direct",
        help=(
            "Use verified ENCODE BEDPE files directly, or filter the precomputed "
            f"association cache (default: direct)."
        ),
    )
    parser.add_argument(
        "--direct-bedpe",
        type=Path,
        nargs="*",
        default=None,
        help=(
            "Released ENCODE BEDPE loop files to use in direct mode "
            f"(default: {' '.join(str(path) for path in DEFAULT_DIRECT_BEDPE)})"
        ),
    )
    parser.add_argument(
        "--ld-source",
        type=Path,
        default=DEFAULT_LD_SOURCE,
        help=f"BED-like SNP source used to place SNP markers (default: {DEFAULT_LD_SOURCE})",
    )
    parser.add_argument(
        "--tag-source",
        type=Path,
        default=DEFAULT_TAG_SOURCE,
        help=f"BED-like tag-SNP source used as an additional coordinate source (default: {DEFAULT_TAG_SOURCE})",
    )
    parser.add_argument(
        "--resolved-snps",
        type=Path,
        default=DEFAULT_RESOLVED_SNPS,
        help=f"BED-like manually/externally resolved SNP coordinates (default: {DEFAULT_RESOLVED_SNPS})",
    )
    parser.add_argument(
        "--provenance",
        type=Path,
        default=DEFAULT_PROVENANCE,
        help=f"JSON provenance metadata for source validation and reporting (default: {DEFAULT_PROVENANCE})",
    )
    parser.add_argument(
        "--output-prefix",
        default=DEFAULT_OUTPUT_PREFIX,
        help=f"Prefix for generated files (default: {DEFAULT_OUTPUT_PREFIX})",
    )
    return parser.parse_args()


def ensure_chr_prefix(chrom: str) -> str:
    value = str(chrom).strip()
    if not value.lower().startswith("chr"):
        value = f"chr{value}"
    return value


def parse_int(value: str) -> int:
    return int(str(value).strip().replace(",", ""))


def split_snps(value: str) -> tuple[str, ...]:
    return tuple(part for part in re.split(r"[\s,;]+", value.strip()) if part)


def snp_key(value: str) -> str:
    return value.strip().lower()


def sanitize_name(value: str, fallback: str = "item", limit: int = 120) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.:|-]+", "_", value.strip())
    cleaned = cleaned.strip("_") or fallback
    return cleaned[:limit]


def evidence_color(evidence: str) -> str:
    lowered = evidence.lower()
    if "kallmann" in lowered:
        return "204,37,41"
    if "paintor" in lowered:
        return "117,112,179"
    return "27,158,119"


def read_regions(path: Path) -> list[Region]:
    regions: list[Region] = []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"region_id", "chrom", "start_1based", "end_1based", "snps"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")

        for line_number, row in enumerate(reader, start=2):
            region_id = (row.get("region_id") or "").strip()
            chrom = ensure_chr_prefix(row.get("chrom") or "")
            snps = split_snps(row.get("snps") or "")
            if not region_id:
                raise ValueError(f"{path}:{line_number}: empty region_id")
            if not snps:
                raise ValueError(f"{path}:{line_number}: empty snps")

            start = parse_int(row.get("start_1based") or "")
            end = parse_int(row.get("end_1based") or "")
            if start <= 0 or end <= 0:
                raise ValueError(f"{path}:{line_number}: coordinates must be positive")
            if end < start:
                start, end = end, start

            regions.append(
                Region(
                    region_id=region_id,
                    chrom=chrom,
                    start_1based=start,
                    end_1based=end,
                    snps=snps,
                    evidence=(row.get("evidence") or "").strip(),
                    note=(row.get("note") or "").strip(),
                )
            )

    if not regions:
        raise ValueError(f"No regions found in {path}")
    return regions


def read_snp_coordinates(path: Path) -> dict[str, SnpCoordinate]:
    coordinates: dict[str, SnpCoordinate] = {}
    if not path.exists():
        return coordinates

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            if not raw_line.startswith("chr"):
                continue
            parts = raw_line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                start0 = int(parts[1])
                end0 = int(parts[2])
            except ValueError:
                continue
            if end0 <= start0:
                end0 = start0 + 1
            coordinates[snp_key(parts[3])] = SnpCoordinate(
                chrom=ensure_chr_prefix(parts[0]),
                start0=start0,
                end0=end0,
            )
    return coordinates


def read_snp_coordinates_from_sources(paths: Sequence[Path]) -> dict[str, SnpCoordinate]:
    coordinates: dict[str, SnpCoordinate] = {}
    for path in paths:
        coordinates.update(read_snp_coordinates(path))
    return coordinates


def parse_region_text(value: str) -> tuple[str, int, int]:
    match = re.fullmatch(r"(chr[0-9A-Za-z]+):([0-9,]+)-([0-9,]+)", value.strip())
    if not match:
        raise ValueError(f"Invalid region: {value}")
    chrom = ensure_chr_prefix(match.group(1))
    start = parse_int(match.group(2))
    end = parse_int(match.group(3))
    if end < start:
        start, end = end, start
    return chrom, start - 1, end


def regions_by_snp(regions: Sequence[Region]) -> dict[str, list[Region]]:
    mapping: dict[str, list[Region]] = {}
    for region in regions:
        for snp in region.snps:
            mapping.setdefault(snp_key(snp), []).append(region)
    return mapping


def intervals_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    return start_a < end_b and start_b < end_a


def read_h1_links(path: Path, regions: Sequence[Region]) -> list[H1Link]:
    if not path.exists():
        raise FileNotFoundError(f"Missing H1-hESC association file: {path}")

    snp_to_regions = regions_by_snp(regions)
    links: list[H1Link] = []
    seen: set[tuple[str, str, str, int, int, str, str]] = set()

    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            query_variant = (row.get("query_variant") or row.get("source_name") or "").strip()
            if not query_variant:
                continue
            matching_regions = snp_to_regions.get(snp_key(query_variant), [])
            if not matching_regions:
                continue

            target_region = (row.get("associated_region") or "").strip()
            if not target_region:
                continue
            try:
                target_chrom, target_start0, target_end0 = parse_region_text(target_region)
            except ValueError:
                continue

            target_name = (row.get("target_name") or row.get("method") or "H1_hESC_link").strip()
            method = (row.get("method") or "").strip()
            for region in matching_regions:
                key = (
                    region.region_id,
                    snp_key(query_variant),
                    target_chrom,
                    target_start0,
                    target_end0,
                    target_name,
                    method,
                )
                if key in seen:
                    continue
                seen.add(key)
                links.append(
                    H1Link(
                        region=region,
                        query_variant=query_variant,
                        target_chrom=target_chrom,
                        target_start0=target_start0,
                        target_end0=target_end0,
                        target_region=target_region,
                        target_name=target_name,
                        method=method,
                    )
                )

    links.sort(
        key=lambda link: (
            link.region.chrom,
            link.region.start0,
            link.query_variant,
            link.target_chrom,
            link.target_start0,
            link.target_end0,
            link.target_name,
        )
    )
    return links


def text_open(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def accession_from_bedpe_path(path: Path) -> str:
    name = path.name
    for suffix in (".bedpe.gz", ".bedpe", ".gz"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def read_direct_bedpe_links(
    paths: Sequence[Path],
    regions: Sequence[Region],
    snp_coordinates: dict[str, SnpCoordinate],
) -> list[H1Link]:
    snp_entries_by_chrom: dict[str, list[tuple[int, int, Region, str]]] = {}
    for region in regions:
        for snp in region.snps:
            coordinate = snp_coordinates.get(snp_key(snp))
            if not coordinate:
                continue
            snp_entries_by_chrom.setdefault(coordinate.chrom, []).append(
                (coordinate.start0, coordinate.end0, region, snp)
            )

    links: list[H1Link] = []
    seen: set[tuple[str, str, str, int, int, str, str]] = set()

    for path in paths:
        if not path.exists():
            raise FileNotFoundError(f"Missing direct BEDPE source: {path}")
        accession = accession_from_bedpe_path(path)
        with text_open(path) as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                if not raw_line.strip() or raw_line.startswith(("#", "track", "browser")):
                    continue
                parts = raw_line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                try:
                    chrom1 = ensure_chr_prefix(parts[0])
                    start1 = int(parts[1])
                    end1 = int(parts[2])
                    chrom2 = ensure_chr_prefix(parts[3])
                    start2 = int(parts[4])
                    end2 = int(parts[5])
                except ValueError:
                    continue
                if end1 <= start1 or end2 <= start2:
                    continue

                target_name = parts[6].strip() if len(parts) > 6 and parts[6].strip() else accession
                method = f"ENCODE:{accession}:line_{line_number}"

                candidates = [
                    (chrom1, start1, end1, chrom2, start2, end2),
                    (chrom2, start2, end2, chrom1, start1, end1),
                ]
                for source_chrom, source_start, source_end, target_chrom, target_start, target_end in candidates:
                    for snp_start, snp_end, region, snp in snp_entries_by_chrom.get(source_chrom, []):
                        if not intervals_overlap(snp_start, snp_end, source_start, source_end):
                            continue
                        key = (
                            region.region_id,
                            snp_key(snp),
                            target_chrom,
                            target_start,
                            target_end,
                            target_name,
                            method,
                        )
                        if key in seen:
                            continue
                        seen.add(key)
                        links.append(
                            H1Link(
                                region=region,
                                query_variant=snp,
                                target_chrom=target_chrom,
                                target_start0=target_start,
                                target_end0=target_end,
                                target_region=format_region_1based(target_chrom, target_start, target_end),
                                target_name=target_name,
                                method=method,
                                source_accession=accession,
                            )
                        )

    links.sort(
        key=lambda link: (
            link.region.chrom,
            link.region.start0,
            link.query_variant,
            link.target_chrom,
            link.target_start0,
            link.target_end0,
            link.target_name,
            link.method,
        )
    )
    return links


def format_region_1based(chrom: str, start0: int, end0: int) -> str:
    return f"{ensure_chr_prefix(chrom)}:{start0 + 1}-{end0}"


def track_header(track_type: str, name: str, description: str, extra: str = "") -> str:
    suffix = f" {extra}" if extra else ""
    return f'track type={track_type} name="{name}" description="{description}"{suffix}'


def bed_track_header(name: str, description: str, extra: str = "") -> str:
    suffix = f" {extra}" if extra else ""
    return f'track name="{name}" description="{description}"{suffix}'


def region_track_lines(regions: Sequence[Region]) -> list[str]:
    lines = [
        bed_track_header(
            "Luciferase_H1_hESC_regions",
            "Luciferase active and PAINTOR regions for H1-hESC review",
            "visibility=pack itemRgb=On",
        )
    ]
    for region in regions:
        snp_label = ",".join(region.snps)
        label = sanitize_name(f"{region.region_id}|{snp_label}", region.region_id)
        lines.append(
            "\t".join(
                [
                    region.chrom,
                    str(region.start0),
                    str(region.end0),
                    label,
                    "1000",
                    ".",
                    str(region.start0),
                    str(region.end0),
                    evidence_color(region.evidence),
                ]
            )
        )
    return lines


def snp_track_lines(regions: Sequence[Region], snp_coordinates: dict[str, SnpCoordinate], links: Sequence[H1Link]) -> list[str]:
    lines = [
        bed_track_header(
            "Luciferase_H1_hESC_SNPs",
            "SNPs listed for luciferase active and PAINTOR regions",
            "visibility=pack itemRgb=On",
        )
    ]
    snps_with_links = {snp_key(link.query_variant) for link in links}
    for region in regions:
        for snp in region.snps:
            coord = snp_coordinates.get(snp_key(snp))
            if not coord:
                continue
            score = "1000" if snp_key(snp) in snps_with_links else "500"
            label = sanitize_name(f"{snp}|{region.region_id}", snp)
            lines.append(
                "\t".join(
                    [
                        coord.chrom,
                        str(coord.start0),
                        str(coord.end0),
                        label,
                        score,
                        ".",
                        str(coord.start0),
                        str(coord.end0),
                        evidence_color(region.evidence),
                    ]
                )
            )
    return lines


def interaction_value(target_name: str) -> str:
    value = target_name.strip()
    if re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", value):
        return value
    return "1"


def interact_track_lines(links: Sequence[H1Link]) -> list[str]:
    lines = [
        track_header(
            "interact",
            "Luciferase_H1_hESC_links",
            "Verified ENCODE H1 ChIA-PET loop arcs for requested SNPs",
            "visibility=full maxHeightPixels=256:128:64",
        )
    ]
    for idx, link in enumerate(links, start=1):
        span_chrom = link.region.chrom
        span_start = link.region.start0
        span_end = link.region.end0
        if link.region.chrom == link.target_chrom:
            span_start = min(link.region.start0, link.target_start0)
            span_end = max(link.region.end0, link.target_end0)

        source_name = sanitize_name(f"{link.region.region_id}|{link.query_variant}", link.region.region_id)
        row_name = sanitize_name(f"{source_name}|H1_{idx}", f"H1_{idx}")
        target_name = sanitize_name(link.target_name, "H1_target")
        lines.append(
            "\t".join(
                [
                    span_chrom,
                    str(span_start),
                    str(span_end),
                    row_name,
                    "1000",
                    interaction_value(link.target_name),
                    ".",
                    evidence_color(link.region.evidence),
                    link.region.chrom,
                    str(link.region.start0),
                    str(link.region.end0),
                    source_name,
                    ".",
                    link.target_chrom,
                    str(link.target_start0),
                    str(link.target_end0),
                    target_name,
                    ".",
                ]
            )
        )
    return lines


def write_lines(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_summary(
    path: Path,
    regions: Sequence[Region],
    snp_coordinates: dict[str, SnpCoordinate],
    links: Sequence[H1Link],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    links_by_region: dict[str, list[H1Link]] = {region.region_id: [] for region in regions}
    for link in links:
        links_by_region.setdefault(link.region.region_id, []).append(link)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "region_id",
                "region",
                "snps",
                "evidence",
                "note",
                "coordinate_snps_found",
                "coordinate_snps_missing",
                "snps_with_h1_interactions",
                "h1_interaction_count",
            ]
        )
        for region in regions:
            region_links = links_by_region.get(region.region_id, [])
            found = [snp for snp in region.snps if snp_key(snp) in snp_coordinates]
            missing = [snp for snp in region.snps if snp_key(snp) not in snp_coordinates]
            snps_with_interactions = sorted({link.query_variant for link in region_links})
            writer.writerow(
                [
                    region.region_id,
                    region.region_text,
                    " ".join(region.snps),
                    region.evidence,
                    region.note,
                    " ".join(found),
                    " ".join(missing),
                    " ".join(snps_with_interactions),
                    len(region_links),
                ]
            )


def write_links_table(path: Path, links: Sequence[H1Link]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "region_id",
                "source_region",
                "query_variant",
                "associated_region",
                "target_name",
                "method",
                "source_accession",
            ]
        )
        for link in links:
            writer.writerow(
                [
                    link.region.region_id,
                    link.region.region_text,
                    link.query_variant,
                    link.target_region,
                    link.target_name,
                    link.method,
                    link.source_accession,
                ]
            )


def write_report(
    path: Path,
    regions: Sequence[Region],
    snp_coordinates: dict[str, SnpCoordinate],
    links: Sequence[H1Link],
    output_paths: dict[str, Path],
    h1_association_path: Path,
    coordinate_sources: Sequence[Path],
    source_mode: str,
    direct_bedpe_paths: Sequence[Path],
    provenance: dict,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    all_snps = [snp for region in regions for snp in region.snps]
    missing_snps = [snp for snp in all_snps if snp_key(snp) not in snp_coordinates]
    snps_with_links = sorted({link.query_variant for link in links})

    links_by_region: dict[str, list[H1Link]] = {region.region_id: [] for region in regions}
    for link in links:
        links_by_region.setdefault(link.region.region_id, []).append(link)

    lines = [
        "# Luciferase H1-hESC Browser Track Report",
        "",
        "## Summary",
        "",
        f"- `regions`: {len(regions)}",
        f"- `snps_requested`: {len(all_snps)}",
        f"- `snps_found_in_coordinate_sources`: {len(all_snps) - len(missing_snps)}",
        f"- `snps_missing_from_coordinate_sources`: {len(missing_snps)}",
        f"- `snps_with_h1_interactions`: {len(snps_with_links)}",
        f"- `h1_interactions`: {len(links)}",
        f"- `source_mode`: {source_mode}",
        f"- `h1_association_source`: {h1_association_path if source_mode == 'cached' else 'direct ENCODE BEDPE files'}",
        "",
        "## Outputs",
        "",
    ]
    for label, output_path in output_paths.items():
        lines.append(f"- `{label}`: `{output_path}`")

    lines.extend(["", "## Source Provenance", ""])
    if provenance:
        lines.extend(
            [
                f"- `metadata_checked_at`: {provenance.get('official_metadata_checked_at', 'unknown')}",
                f"- `metadata_source`: {provenance.get('official_metadata_source', 'unknown')}",
                f"- `interaction_source`: {provenance.get('interaction_source_assertion', 'unknown')}",
                f"- `cell_line_scope`: {provenance.get('cell_line_scope', 'unknown')}",
                f"- `assay_scope`: {provenance.get('assay_scope', 'unknown')}",
                f"- `genome_assembly`: {provenance.get('genome_assembly', 'unknown')}",
            ]
        )
        encode_files = provenance.get("encode_files") or []
        if encode_files:
            lines.extend(
                [
                    "",
                    "| accession | experiment | assay/target | status | assembly | md5 |",
                    "|---|---|---|---|---|---|",
                ]
            )
            for item in encode_files:
                lines.append(
                    "| "
                    + " | ".join(
                        [
                            str(item.get("accession", "")),
                            str(item.get("experiment", "")),
                            f"{item.get('experiment_description', '')}; {item.get('target', '')}",
                            str(item.get("status", "")),
                            str(item.get("assembly", "")),
                            str(item.get("md5sum", "")),
                        ]
                    )
                    + " |"
                )
    else:
        lines.append("- No provenance JSON was loaded.")

    lines.extend(["", "## Input Coordinate Sources", ""])
    for coordinate_source in coordinate_sources:
        lines.append(f"- `{coordinate_source}`")
    if direct_bedpe_paths:
        lines.extend(["", "## Direct BEDPE Sources", ""])
        for direct_bedpe_path in direct_bedpe_paths:
            lines.append(f"- `{direct_bedpe_path}`")

    lines.extend(
        [
            "",
            "## Per-Region H1-hESC Interactions",
            "",
            "| region_id | region | snps | evidence | h1_interactions | snps_with_interactions | note |",
            "|---|---|---|---|---:|---|---|",
        ]
    )
    for region in regions:
        region_links = links_by_region.get(region.region_id, [])
        snps_for_region = sorted({link.query_variant for link in region_links})
        lines.append(
            "| "
            + " | ".join(
                [
                    region.region_id,
                    region.region_text,
                    " ".join(region.snps),
                    region.evidence,
                    str(len(region_links)),
                    " ".join(snps_for_region) if snps_for_region else "none",
                    region.note,
                ]
            )
            + " |"
        )

    lines.extend(["", "## Missing From Coordinate Sources", ""])
    if missing_snps:
        lines.append("- " + ", ".join(missing_snps))
    else:
        lines.append("- None")
    lines.extend(
        [
            "",
            "## UCSC Use",
            "",
            f"Upload `{output_paths['ucsc_session']}` as a custom track on hg38. "
            "The combined file includes region, SNP, and H1-hESC interaction tracks.",
        ]
    )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_provenance(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    args = parse_args()
    regions = read_regions(args.regions)
    coordinate_sources = [args.ld_source, args.tag_source, args.resolved_snps]
    snp_coordinates = read_snp_coordinates_from_sources(coordinate_sources)
    direct_bedpe_paths = list(args.direct_bedpe) if args.direct_bedpe is not None else list(DEFAULT_DIRECT_BEDPE)
    if args.source_mode == "direct":
        links = read_direct_bedpe_links(direct_bedpe_paths, regions, snp_coordinates)
    else:
        links = read_h1_links(args.h1_associations, regions)
    provenance = read_provenance(args.provenance)

    prefix = Path(args.output_prefix)
    output_paths = {
        "regions_bed": prefix.with_name(prefix.name + "_regions.bed"),
        "snps_bed": prefix.with_name(prefix.name + "_snps.bed"),
        "links_interact": prefix.with_name(prefix.name + "_links.interact"),
        "links_tsv": prefix.with_name(prefix.name + "_links.tsv"),
        "ucsc_session": prefix.with_name(prefix.name + "_ucsc_session.txt"),
        "summary_tsv": prefix.with_name(prefix.name + "_summary.tsv"),
        "report_md": prefix.with_name(prefix.name + "_report.md"),
    }

    region_lines = region_track_lines(regions)
    snp_lines = snp_track_lines(regions, snp_coordinates, links)
    interact_lines = interact_track_lines(links)

    write_lines(output_paths["regions_bed"], region_lines)
    write_lines(output_paths["snps_bed"], snp_lines)
    write_lines(output_paths["links_interact"], interact_lines)
    write_links_table(output_paths["links_tsv"], links)

    first_region = regions[0]
    session_lines = [
        f"browser position {first_region.region_text}",
        *region_lines,
        *snp_lines,
        *interact_lines,
    ]
    write_lines(output_paths["ucsc_session"], session_lines)
    write_summary(output_paths["summary_tsv"], regions, snp_coordinates, links)
    write_report(
        output_paths["report_md"],
        regions,
        snp_coordinates,
        links,
        output_paths,
        args.h1_associations,
        coordinate_sources,
        args.source_mode,
        direct_bedpe_paths if args.source_mode == "direct" else [],
        provenance,
    )

    all_snps = [snp for region in regions for snp in region.snps]
    missing_snps = [snp for snp in all_snps if snp_key(snp) not in snp_coordinates]
    print(
        "Wrote focused H1-hESC browser tracks "
        f"for {len(regions)} regions, {len(all_snps)} SNPs, and {len(links)} interactions "
        f"({len(missing_snps)} SNPs missing from coordinate sources; source_mode={args.source_mode})."
    )
    for label, output_path in output_paths.items():
        print(f"{label}: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
