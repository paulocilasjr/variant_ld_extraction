#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import time
from collections import defaultdict
from pathlib import Path

import requests

GWAS_REGION_URL = "https://www.ebi.ac.uk/gwas/api/v2/regions/{region}/associations"
LDPAIR_URL = "https://ldlink.nih.gov/LDlinkRestWeb/ldpair"
DEFAULT_INPUT = Path("TECAC_GWAS_index_SNPs_OCT_2025")
DEFAULT_OUTPUT_PREFIX = "tag_snp_2mb"
DEFAULT_WINDOW_BP = 1_000_000
PAGE_SIZE = 500

ALL_FIELDS = [
    "tag_snp",
    "tag_chrom",
    "tag_pos",
    "source_region",
    "window_start",
    "window_end",
    "association_variant",
    "association_label",
    "association_chrom",
    "association_pos",
    "distance_to_tag_bp",
    "trait_name",
    "trait_category",
    "matched_keywords",
    "mapped_genes",
    "accession_id",
    "pubmed_id",
    "author",
    "publication_date",
    "p_value",
    "beta",
    "ci",
    "risk_frequency",
]

RELEVANT_FIELDS = ALL_FIELDS + [
    "ldpair_r2",
    "ldpair_d_prime",
    "ldpair_p",
    "ldpair_comment",
]

SUMMARY_FIELDS = [
    "tag_snp",
    "tag_chrom",
    "tag_pos",
    "source_region",
    "window_start",
    "window_end",
    "all_rows",
    "relevant_rows",
    "testicular_cancer_rows",
    "other_cancer_rows",
    "testicular_related_rows",
    "top_testicular_cancer_r2",
    "top_testicular_cancer_trait",
    "top_testicular_cancer_variant",
    "top_other_cancer_r2",
    "top_other_cancer_trait",
    "top_other_cancer_variant",
    "top_testicular_related_r2",
    "top_testicular_related_trait",
    "top_testicular_related_variant",
]

TESTICULAR_TERMS = {
    "testicular",
    "testis",
    "germ cell tumor",
    "germ cell tumour",
    "seminoma",
    "nonseminoma",
    "cryptorchid",
    "hypospad",
    "azoosperm",
    "oligosperm",
    "sperm",
    "semen",
    "infertility",
    "fertility",
    "testosterone",
    "androgen",
    "gonadal",
    "puberty",
    "voice breaking",
    "facial hair",
}

CANCER_TERMS = {
    "cancer",
    "carcinoma",
    "tumor",
    "tumour",
    "neoplasm",
    "leukemia",
    "lymphoma",
    "melanoma",
    "sarcoma",
    "myeloma",
    "glioma",
    "adenoma",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Scan a 2 Mb tag-SNP-centered GWAS window, retain cancer and "
            "testicular-related traits, and quantify pairwise LD to each tag SNP."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"BED-like tag SNP file (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--output-prefix",
        default=DEFAULT_OUTPUT_PREFIX,
        help=f"Output prefix (default: {DEFAULT_OUTPUT_PREFIX})",
    )
    parser.add_argument(
        "--window-bp",
        type=int,
        default=DEFAULT_WINDOW_BP,
        help="Number of base pairs on each side of the tag SNP (default: 1000000).",
    )
    parser.add_argument(
        "--regions-file",
        type=Path,
        help="Optional BED-style regions file; only tag SNPs overlapping these regions are scanned.",
    )
    parser.add_argument(
        "--sleep-seconds",
        type=float,
        default=0.15,
        help="Delay between LDpair requests (default: 0.15).",
    )
    parser.add_argument(
        "--max-tag-snps",
        type=int,
        default=0,
        help="Optional cap for test runs; 0 means all tag SNPs.",
    )
    return parser.parse_args()


def scientific_pvalue(base: object, exponent: object) -> str:
    if base in (None, ""):
        return ""
    base_str = str(base).strip()
    exponent_str = str(exponent).strip()
    if exponent_str:
        return f"{base_str}e{exponent_str}"
    return base_str


def normalize_trait_name(association: dict) -> str:
    trait_name = association.get("traitName")
    if isinstance(trait_name, list):
        values = [str(value).strip() for value in trait_name if str(value).strip()]
        if values:
            return "; ".join(dict.fromkeys(values))
    elif isinstance(trait_name, str) and trait_name.strip():
        return trait_name.strip()

    labels = []
    for trait in association.get("efoTraits", []) or []:
        if isinstance(trait, dict):
            label = str(trait.get("label", "")).strip()
            if label:
                labels.append(label)
    return "; ".join(dict.fromkeys(labels))


def extract_risk_alleles(association: dict) -> list[dict]:
    values = association.get("riskAllele")
    if values is None:
        return []
    candidates = values if isinstance(values, list) else [values]

    rows = []
    seen = set()
    for candidate in candidates:
        if isinstance(candidate, dict):
            key = str(candidate.get("key", "")).strip()
            label = str(candidate.get("label", "")).strip() or key
        else:
            key = str(candidate).strip()
            label = key
        if not key.startswith("rs"):
            continue
        marker = (key, label)
        if marker in seen:
            continue
        seen.add(marker)
        rows.append({"key": key, "label": label})
    return rows


def parse_location(value: str) -> tuple[str, int] | tuple[None, None]:
    if not isinstance(value, str) or ":" not in value:
        return None, None
    chrom, pos_str = value.split(":", 1)
    try:
        return f"chr{chrom.removeprefix('chr')}", int(pos_str.replace(",", ""))
    except ValueError:
        return None, None


def classify_trait(trait_name: str) -> tuple[str, list[str]]:
    lowered = trait_name.lower()
    matched_testicular = sorted(term for term in TESTICULAR_TERMS if term in lowered)
    matched_cancer = sorted(term for term in CANCER_TERMS if term in lowered)

    if matched_testicular and matched_cancer:
        return "testicular_cancer", matched_testicular + matched_cancer
    if matched_cancer:
        return "other_cancer", matched_cancer
    if matched_testicular:
        return "testicular_related_phenotype", matched_testicular
    return "other", []


def load_tag_snps(path: Path) -> list[dict]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("chr"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0].strip()
            rsid = parts[3].strip()
            if not rsid.startswith("rs"):
                continue
            try:
                start0 = int(parts[1])
            except ValueError:
                continue
            rows.append(
                {
                    "tag_snp": rsid,
                    "chrom": chrom,
                    "pos": start0 + 1,
                    "source_region": "",
                }
            )
    return rows


def load_regions(path: Path) -> list[dict]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("chr"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                start0 = int(parts[1])
                end0 = int(parts[2])
            except ValueError:
                continue
            if end0 <= start0:
                continue
            rows.append(
                {
                    "chrom": parts[0].strip(),
                    "start0": start0,
                    "end0": end0,
                    "source_region": f"{parts[0].strip()}:{start0 + 1}-{end0}",
                }
            )
    return rows


def restrict_tag_snps_to_regions(tag_snps: list[dict], regions: list[dict]) -> list[dict]:
    filtered = []
    for tag in tag_snps:
        pos0 = int(tag["pos"]) - 1
        for region in regions:
            if (
                tag["chrom"] == region["chrom"]
                and region["start0"] <= pos0 < region["end0"]
            ):
                filtered.append({**tag, "source_region": region["source_region"]})
                break
    return filtered


def fetch_region_associations(
    session: requests.Session,
    chrom: str,
    pos: int,
    window_bp: int,
) -> tuple[list[dict], int, int]:
    start = max(1, pos - window_bp)
    end = pos + window_bp
    api_region = f"{chrom.removeprefix('chr')}:{start}-{end}"
    rows: list[dict] = []
    page = 0

    while True:
        response = session.get(
            GWAS_REGION_URL.format(region=api_region),
            params={"page": page, "size": PAGE_SIZE},
            timeout=60,
        )
        response.raise_for_status()
        payload = response.json()
        current = payload.get("_embedded", {}).get("associations", [])
        if not current:
            break
        rows.extend(current)

        total_pages = payload.get("page", {}).get("totalPages")
        if isinstance(total_pages, int):
            if page + 1 >= total_pages:
                break
        elif len(current) < PAGE_SIZE:
            break
        page += 1

    return rows, start, end


def parse_ldpair_payload(payload: object, tag_snp: str, other_snp: str) -> dict[str, str]:
    if tag_snp == other_snp:
        return {
            "ldpair_r2": "1",
            "ldpair_d_prime": "1",
            "ldpair_p": "",
            "ldpair_comment": "Same SNP as the tag SNP",
        }

    if not isinstance(payload, list) or not payload or not isinstance(payload[0], dict):
        return {
            "ldpair_r2": "",
            "ldpair_d_prime": "",
            "ldpair_p": "",
            "ldpair_comment": f"Unexpected LDpair payload type: {type(payload).__name__}",
        }

    first = payload[0]
    stats = first.get("statistics", {}) if isinstance(first.get("statistics"), dict) else {}
    comments = first.get("corr_alleles", [])
    if isinstance(comments, list):
        comment_text = "; ".join(str(value).strip() for value in comments if str(value).strip())
    else:
        comment_text = str(comments).strip()
    return {
        "ldpair_r2": str(stats.get("r2", "")).strip(),
        "ldpair_d_prime": str(stats.get("d_prime", "")).strip(),
        "ldpair_p": str(stats.get("p", "")).strip(),
        "ldpair_comment": comment_text,
    }


def fetch_ldpair(
    session: requests.Session,
    tag_snp: str,
    other_snp: str,
    sleep_seconds: float,
    cache: dict[tuple[str, str], dict[str, str]],
) -> dict[str, str]:
    key = tuple(sorted((tag_snp, other_snp)))
    if key in cache:
        return cache[key]

    response = session.get(
        LDPAIR_URL,
        params={
            "var1": tag_snp,
            "var2": other_snp,
            "pop": "EUR",
            "genome_build": "grch38",
        },
        timeout=120,
    )
    response.raise_for_status()
    payload = response.json()
    parsed = parse_ldpair_payload(payload, tag_snp, other_snp)
    cache[key] = parsed
    if sleep_seconds > 0:
        time.sleep(sleep_seconds)
    return parsed


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        if rows:
            writer.writerows(rows)


def float_or_negative_inf(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("-inf")


def summarize_tag_rows(tag_rows: list[dict]) -> dict[str, str]:
    summary = {
        "tag_snp": "",
        "tag_chrom": "",
        "tag_pos": "",
        "source_region": "",
        "window_start": "",
        "window_end": "",
        "all_rows": "0",
        "relevant_rows": "0",
        "testicular_cancer_rows": "0",
        "other_cancer_rows": "0",
        "testicular_related_rows": "0",
        "top_testicular_cancer_r2": "",
        "top_testicular_cancer_trait": "",
        "top_testicular_cancer_variant": "",
        "top_other_cancer_r2": "",
        "top_other_cancer_trait": "",
        "top_other_cancer_variant": "",
        "top_testicular_related_r2": "",
        "top_testicular_related_trait": "",
        "top_testicular_related_variant": "",
    }
    if not tag_rows:
        return summary

    first = tag_rows[0]
    summary.update(
        {
            "tag_snp": first["tag_snp"],
            "tag_chrom": first["tag_chrom"],
            "tag_pos": str(first["tag_pos"]),
            "source_region": first.get("source_region", ""),
            "window_start": str(first["window_start"]),
            "window_end": str(first["window_end"]),
            "all_rows": str(len(tag_rows)),
        }
    )

    relevant = [row for row in tag_rows if row["trait_category"] != "other"]
    summary["relevant_rows"] = str(len(relevant))
    summary["testicular_cancer_rows"] = str(
        sum(1 for row in relevant if row["trait_category"] == "testicular_cancer")
    )
    summary["other_cancer_rows"] = str(
        sum(1 for row in relevant if row["trait_category"] == "other_cancer")
    )
    summary["testicular_related_rows"] = str(
        sum(1 for row in relevant if row["trait_category"] == "testicular_related_phenotype")
    )

    for category, prefix in (
        ("testicular_cancer", "top_testicular_cancer"),
        ("other_cancer", "top_other_cancer"),
        ("testicular_related_phenotype", "top_testicular_related"),
    ):
        category_rows = [row for row in relevant if row["trait_category"] == category]
        if not category_rows:
            continue
        best = max(category_rows, key=lambda row: float_or_negative_inf(row.get("ldpair_r2", "")))
        summary[f"{prefix}_r2"] = best.get("ldpair_r2", "")
        summary[f"{prefix}_trait"] = best["trait_name"]
        summary[f"{prefix}_variant"] = best["association_variant"]

    return summary


def main() -> int:
    args = parse_args()
    if not args.input.exists():
        raise FileNotFoundError(f"Missing input file: {args.input}")

    tag_snps = load_tag_snps(args.input)
    if not tag_snps:
        raise RuntimeError(f"No tag SNP rows were loaded from {args.input}")
    if args.regions_file is not None:
        if not args.regions_file.exists():
            raise FileNotFoundError(f"Missing regions file: {args.regions_file}")
        tag_snps = restrict_tag_snps_to_regions(tag_snps, load_regions(args.regions_file))
        if not tag_snps:
            raise RuntimeError(
                f"No tag SNPs from {args.input} overlapped the regions in {args.regions_file}"
            )
    if args.max_tag_snps > 0:
        tag_snps = tag_snps[: args.max_tag_snps]

    output_prefix = Path(args.output_prefix)
    all_rows: list[dict] = []
    relevant_rows: list[dict] = []
    summary_rows: list[dict] = []
    ldpair_cache: dict[tuple[str, str], dict[str, str]] = {}

    with requests.Session() as session:
        for index, tag in enumerate(tag_snps, start=1):
            tag_rows: list[dict] = []
            associations, window_start, window_end = fetch_region_associations(
                session,
                tag["chrom"],
                tag["pos"],
                args.window_bp,
            )
            print(
                f"[{index}/{len(tag_snps)}] {tag['tag_snp']} "
                f"{tag['chrom']}:{window_start}-{window_end} -> {len(associations)} associations"
            )

            seen_rows = set()
            for association in associations:
                trait_name = normalize_trait_name(association)
                trait_category, matched_keywords = classify_trait(trait_name)
                mapped_genes = "; ".join(
                    dict.fromkeys(
                        str(gene).strip()
                        for gene in association.get("mappedGenes", []) or []
                        if str(gene).strip()
                    )
                )

                locations = association.get("locations", []) or []
                association_chrom, association_pos = None, None
                for location in locations:
                    association_chrom, association_pos = parse_location(location)
                    if association_chrom and association_pos:
                        break

                for risk_allele in extract_risk_alleles(association):
                    row = {
                        "tag_snp": tag["tag_snp"],
                        "tag_chrom": tag["chrom"],
                        "tag_pos": tag["pos"],
                        "source_region": tag.get("source_region", ""),
                        "window_start": window_start,
                        "window_end": window_end,
                        "association_variant": risk_allele["key"],
                        "association_label": risk_allele["label"],
                        "association_chrom": association_chrom or "",
                        "association_pos": association_pos or "",
                        "distance_to_tag_bp": (
                            abs(association_pos - tag["pos"]) if association_pos else ""
                        ),
                        "trait_name": trait_name,
                        "trait_category": trait_category,
                        "matched_keywords": "; ".join(matched_keywords),
                        "mapped_genes": mapped_genes,
                        "accession_id": str(association.get("accessionId", "")).strip(),
                        "pubmed_id": str(association.get("pubmedId", "")).strip(),
                        "author": str(association.get("author", "")).strip(),
                        "publication_date": str(association.get("publicationDate", "")).strip(),
                        "p_value": scientific_pvalue(
                            association.get("pValue"),
                            association.get("pValueExponent"),
                        ),
                        "beta": str(association.get("beta", "")).strip(),
                        "ci": str(association.get("ci", "")).strip(),
                        "risk_frequency": str(association.get("riskFrequency", "")).strip(),
                    }
                    marker = tuple(str(row[field]) for field in ALL_FIELDS)
                    if marker in seen_rows:
                        continue
                    seen_rows.add(marker)
                    all_rows.append(row)

                    if trait_category == "other":
                        tag_rows.append(row)
                        continue

                    try:
                        ld_metrics = fetch_ldpair(
                            session,
                            tag["tag_snp"],
                            risk_allele["key"],
                            args.sleep_seconds,
                            ldpair_cache,
                        )
                    except requests.RequestException as exc:
                        ld_metrics = {
                            "ldpair_r2": "",
                            "ldpair_d_prime": "",
                            "ldpair_p": "",
                            "ldpair_comment": f"LDpair request failed: {exc}",
                        }

                    relevant_row = {**row, **ld_metrics}
                    relevant_rows.append(relevant_row)
                    tag_rows.append(relevant_row)

            summary_rows.append(summarize_tag_rows(tag_rows))

    all_rows.sort(key=lambda row: (row["tag_chrom"], int(row["tag_pos"]), row["trait_name"]))
    relevant_rows.sort(
        key=lambda row: (
            row["tag_chrom"],
            int(row["tag_pos"]),
            row["trait_category"],
            -float_or_negative_inf(row.get("ldpair_r2", "")),
            row["trait_name"],
        )
    )
    summary_rows.sort(key=lambda row: (row["tag_chrom"], int(row["tag_pos"] or 0), row["tag_snp"]))

    write_tsv(Path(f"{output_prefix}_all_associations.tsv"), ALL_FIELDS, all_rows)
    write_tsv(Path(f"{output_prefix}_relevant_associations.tsv"), RELEVANT_FIELDS, relevant_rows)
    write_tsv(Path(f"{output_prefix}_summary.tsv"), SUMMARY_FIELDS, summary_rows)

    report = {
        "input": str(args.input),
        "regions_file": str(args.regions_file) if args.regions_file is not None else "",
        "window_bp_each_side": args.window_bp,
        "tag_snp_count": len(tag_snps),
        "all_rows": len(all_rows),
        "relevant_rows": len(relevant_rows),
        "ldpair_cache_entries": len(ldpair_cache),
    }
    Path(f"{output_prefix}_report.json").write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )

    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
