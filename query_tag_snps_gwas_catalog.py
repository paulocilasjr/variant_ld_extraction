#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import requests

GWAS_VARIANT_URL = "https://www.ebi.ac.uk/gwas/api/v2/variants/{rsid}/associations"
DEFAULT_INPUT = Path("TECAC_GWAS_index_SNPs_OCT_2025")
DEFAULT_OUTPUT_PREFIX = "tag_snp_gwas_catalog"
PAGE_SIZE = 500

ALL_FIELDS = [
    "tag_snp",
    "tag_chrom",
    "tag_pos",
    "source_region",
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
    "or_value",
    "ci",
    "risk_frequency",
    "association_locations",
]

SUMMARY_FIELDS = [
    "tag_snp",
    "tag_chrom",
    "tag_pos",
    "source_region",
    "all_association_rows",
    "cancer_rows",
    "testicular_cancer_rows",
    "other_cancer_rows",
    "top_testicular_cancer_trait",
    "top_testicular_cancer_p_value",
    "top_other_cancer_trait",
    "top_other_cancer_p_value",
]

TESTICULAR_TERMS = {
    "testicular",
    "testis",
    "germ cell tumor",
    "germ cell tumour",
    "seminoma",
    "nonseminoma",
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
        description="Query GWAS Catalog exact tag-SNP associations and retain cancer matches."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument(
        "--regions-file",
        type=Path,
        help="Optional BED-style regions file; only tag SNPs overlapping these regions are queried.",
    )
    parser.add_argument("--output-prefix", default=DEFAULT_OUTPUT_PREFIX)
    return parser.parse_args()


def scientific_pvalue(base: object, exponent: object) -> str:
    if base in (None, ""):
        return ""
    base_str = str(base).strip()
    exponent_str = str(exponent).strip()
    return f"{base_str}e{exponent_str}" if exponent_str else base_str


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


def classify_trait(trait_name: str) -> tuple[str, list[str]]:
    lowered = trait_name.lower()
    matched_testicular = sorted(term for term in TESTICULAR_TERMS if term in lowered)
    matched_cancer = sorted(term for term in CANCER_TERMS if term in lowered)
    if matched_testicular and matched_cancer:
        return "testicular_cancer", matched_testicular + matched_cancer
    if matched_cancer:
        return "other_cancer", matched_cancer
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
                    "tag_chrom": chrom,
                    "tag_pos": start0 + 1,
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
        pos0 = int(tag["tag_pos"]) - 1
        for region in regions:
            if (
                tag["tag_chrom"] == region["chrom"]
                and region["start0"] <= pos0 < region["end0"]
            ):
                filtered.append({**tag, "source_region": region["source_region"]})
                break
    return filtered


def fetch_variant_associations(session: requests.Session, rsid: str) -> list[dict]:
    rows: list[dict] = []
    page = 0
    while True:
        response = session.get(
            GWAS_VARIANT_URL.format(rsid=rsid),
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
    return rows


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        if rows:
            writer.writerows(rows)


def parse_scientific_as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("inf")


def summarize_tag_rows(tag: dict, rows: list[dict]) -> dict[str, str]:
    summary = {
        "tag_snp": tag["tag_snp"],
        "tag_chrom": tag["tag_chrom"],
        "tag_pos": str(tag["tag_pos"]),
        "source_region": tag.get("source_region", ""),
        "all_association_rows": str(len(rows)),
        "cancer_rows": str(sum(1 for row in rows if row["trait_category"] != "other")),
        "testicular_cancer_rows": str(
            sum(1 for row in rows if row["trait_category"] == "testicular_cancer")
        ),
        "other_cancer_rows": str(
            sum(1 for row in rows if row["trait_category"] == "other_cancer")
        ),
        "top_testicular_cancer_trait": "",
        "top_testicular_cancer_p_value": "",
        "top_other_cancer_trait": "",
        "top_other_cancer_p_value": "",
    }

    for category, trait_key, p_key in (
        ("testicular_cancer", "top_testicular_cancer_trait", "top_testicular_cancer_p_value"),
        ("other_cancer", "top_other_cancer_trait", "top_other_cancer_p_value"),
    ):
        category_rows = [row for row in rows if row["trait_category"] == category]
        if not category_rows:
            continue
        best = min(category_rows, key=lambda row: parse_scientific_as_float(row["p_value"]))
        summary[trait_key] = best["trait_name"]
        summary[p_key] = best["p_value"]

    return summary


def main() -> int:
    args = parse_args()
    if not args.input.exists():
        raise FileNotFoundError(f"Missing input file: {args.input}")

    tag_snps = load_tag_snps(args.input)
    if args.regions_file is not None:
        if not args.regions_file.exists():
            raise FileNotFoundError(f"Missing regions file: {args.regions_file}")
        tag_snps = restrict_tag_snps_to_regions(tag_snps, load_regions(args.regions_file))
    if not tag_snps:
        raise RuntimeError("No tag SNPs selected for GWAS Catalog querying")

    all_rows: list[dict] = []
    cancer_rows: list[dict] = []
    summary_rows: list[dict] = []

    with requests.Session() as session:
        for idx, tag in enumerate(tag_snps, start=1):
            associations = fetch_variant_associations(session, tag["tag_snp"])
            print(f"[{idx}/{len(tag_snps)}] {tag['tag_snp']} -> {len(associations)} associations")
            tag_rows = []
            for association in associations:
                trait_name = normalize_trait_name(association)
                trait_category, matched_keywords = classify_trait(trait_name)
                row = {
                    "tag_snp": tag["tag_snp"],
                    "tag_chrom": tag["tag_chrom"],
                    "tag_pos": tag["tag_pos"],
                    "source_region": tag.get("source_region", ""),
                    "trait_name": trait_name,
                    "trait_category": trait_category,
                    "matched_keywords": "; ".join(matched_keywords),
                    "mapped_genes": "; ".join(
                        dict.fromkeys(
                            str(gene).strip()
                            for gene in association.get("mappedGenes", []) or []
                            if str(gene).strip()
                        )
                    ),
                    "accession_id": str(association.get("accessionId", "")).strip(),
                    "pubmed_id": str(association.get("pubmedId", "")).strip(),
                    "author": str(association.get("author", "")).strip(),
                    "publication_date": str(association.get("publicationDate", "")).strip(),
                    "p_value": scientific_pvalue(
                        association.get("pValue"),
                        association.get("pValueExponent"),
                    ),
                    "beta": str(association.get("beta", "")).strip(),
                    "or_value": str(association.get("orValue", "")).strip(),
                    "ci": str(association.get("ci", "")).strip(),
                    "risk_frequency": str(association.get("riskFrequency", "")).strip(),
                    "association_locations": "; ".join(
                        str(location).strip()
                        for location in association.get("locations", []) or []
                        if str(location).strip()
                    ),
                }
                all_rows.append(row)
                tag_rows.append(row)
                if trait_category != "other":
                    cancer_rows.append(row)
            summary_rows.append(summarize_tag_rows(tag, tag_rows))

    output_prefix = Path(args.output_prefix)
    summary_rows.sort(key=lambda row: (row["tag_chrom"], int(row["tag_pos"]), row["tag_snp"]))
    cancer_rows.sort(key=lambda row: (row["tag_chrom"], int(row["tag_pos"]), row["trait_name"]))
    all_rows.sort(key=lambda row: (row["tag_chrom"], int(row["tag_pos"]), row["trait_name"]))

    write_tsv(Path(f"{output_prefix}_all.tsv"), ALL_FIELDS, all_rows)
    write_tsv(Path(f"{output_prefix}_cancer.tsv"), ALL_FIELDS, cancer_rows)
    write_tsv(Path(f"{output_prefix}_summary.tsv"), SUMMARY_FIELDS, summary_rows)

    report = {
        "input": str(args.input),
        "regions_file": str(args.regions_file) if args.regions_file else "",
        "tag_snp_count": len(tag_snps),
        "all_rows": len(all_rows),
        "cancer_rows": len(cancer_rows),
    }
    Path(f"{output_prefix}_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
