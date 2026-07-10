#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


DEFAULT_REGIONS = Path("luciferase_hesc_h1_regions.tsv")
DEFAULT_H1_ASSOCIATIONS = Path("ld_r2_equal_higher_0.8_hESC_H1_associated_regions.tsv")
DEFAULT_LD_SOURCE = Path("ld_r2_equal_higher_0.8")
DEFAULT_OUTPUT_PREFIX = "luciferase_hesc_h1"


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create focused UCSC Genome Browser tracks for luciferase-active "
            "regions evaluated against cached H1-hESC chromatin interactions."
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
        "--ld-source",
        type=Path,
        default=DEFAULT_LD_SOURCE,
        help=f"BED-like SNP source used to place SNP markers (default: {DEFAULT_LD_SOURCE})",
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
            "Cached H1-hESC chromatin interactions for requested SNPs",
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
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_summary(
    path: Path,
    regions: Sequence[Region],
    snp_coordinates: dict[str, SnpCoordinate],
    links: Sequence[H1Link],
) -> None:
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
                "cached_h1_snps_found",
                "cached_h1_snps_missing",
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


def write_report(
    path: Path,
    regions: Sequence[Region],
    snp_coordinates: dict[str, SnpCoordinate],
    links: Sequence[H1Link],
    output_paths: dict[str, Path],
    h1_association_path: Path,
) -> None:
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
        f"- `snps_found_in_cached_source`: {len(all_snps) - len(missing_snps)}",
        f"- `snps_missing_from_cached_source`: {len(missing_snps)}",
        f"- `snps_with_cached_h1_interactions`: {len(snps_with_links)}",
        f"- `h1_interactions`: {len(links)}",
        f"- `h1_association_source`: {h1_association_path}",
        "",
        "## Outputs",
        "",
    ]
    for label, output_path in output_paths.items():
        lines.append(f"- `{label}`: `{output_path}`")

    lines.extend(
        [
            "",
            "## Per-Region H1-hESC Interactions",
            "",
            "| region_id | region | snps | evidence | cached_h1_interactions | snps_with_interactions | note |",
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

    lines.extend(["", "## Missing From Cached H1 Source", ""])
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


def main() -> int:
    args = parse_args()
    regions = read_regions(args.regions)
    snp_coordinates = read_snp_coordinates(args.ld_source)
    links = read_h1_links(args.h1_associations, regions)

    prefix = Path(args.output_prefix)
    output_paths = {
        "regions_bed": prefix.with_name(prefix.name + "_regions.bed"),
        "snps_bed": prefix.with_name(prefix.name + "_snps.bed"),
        "links_interact": prefix.with_name(prefix.name + "_links.interact"),
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

    first_region = regions[0]
    session_lines = [
        f"browser position {first_region.region_text}",
        *region_lines,
        *snp_lines,
        *interact_lines,
    ]
    write_lines(output_paths["ucsc_session"], session_lines)
    write_summary(output_paths["summary_tsv"], regions, snp_coordinates, links)
    write_report(output_paths["report_md"], regions, snp_coordinates, links, output_paths, args.h1_associations)

    all_snps = [snp for region in regions for snp in region.snps]
    missing_snps = [snp for snp in all_snps if snp_key(snp) not in snp_coordinates]
    print(
        "Wrote focused H1-hESC browser tracks "
        f"for {len(regions)} regions, {len(all_snps)} SNPs, and {len(links)} interactions "
        f"({len(missing_snps)} SNPs missing from cached source)."
    )
    for label, output_path in output_paths.items():
        print(f"{label}: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
