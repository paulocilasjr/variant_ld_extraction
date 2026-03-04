#!/usr/bin/env python3
import csv
from pathlib import Path

import requests

from query_ld_link import fetch_ldtrait_rows

region_list = [
    "chr4:103,346,459-103,348,483",
    "chr4:103,438,488-103,440,510",
    "chr4:103,450,003-103,452,005",
    "chr16:74,695,852-74,697,881",
    "chr22:41,021,294-41,023,351",
    "chr22:41,035,391-41,037,269",
    "chr22:41,090,314-41,092,353",
    "chr22:41,099,764-41,101,779",
    "chr22:41,169,600-41,171,531",
    "chr22:41,195,746-41,197,756",
    "chr22:41,216,563-41,218,644",
    "chr22:41,221,933-41,223,980",
    "chr22:41,239,974-41,241,988",
    "chr7:21,455,713-21,457,713",
    "chr19:36,256,494-36,258,488",
    "chr19:36,265,454-36,267,498",
    "chr19:36,308,549-36,310,552",
    "chr19:36,324,138-36,326,171",
    "chr19:36,330,894-36,332,896",
    "chr19:36,342,158-36,344,161",
    "chr19:36,377,852-36,379,868",
    "chr19:36,417,767-36,419,761",
    "chr8:94,615,501-94,617,498",
    "chr3:141,910,334-141,912,334",
    "chr3:141,920,302-141,922,312",
    "chr12:31,949,862-31,951,855",
    "chr9:124,560,615-124,562,618",
    "chr9:124,572,929-124,574,941",
    "chr9:124,593,082-124,594,989",
    "chr9:124,601,306-124,603,321",
    "chr9:124,620,553-124,622,642",
    "chr9:124,659,009-124,661,009",
    "chr9:124,695,030-124,697,031",
    "chr9:124,697,464-124,699,474",
    "chr9:124,733,875-124,735,971",
    "chr9:124,757,413-124,759,380",
    "chr9:124,769,804-124,771,893",
    "chr9:124,794,241-124,796,251",
]

ASSOCIATIONS_URL = "https://www.ebi.ac.uk/gwas/api/v2/regions/{region}/associations"
WINDOW_BP = 500000
PAGE_SIZE = 500
OUTPUT_CSV = Path("risk_allele_tag_ld_matches.csv")
OUTPUT_FIELDS = ["riskAllele", "label", "traitName", "tagSNP", "r2"]


def parse_region(region: str) -> tuple[str, int, int]:
    try:
        chrom, coordinates = region.split(":", 1)
        start_str, end_str = coordinates.split("-", 1)
        start = int(start_str.replace(",", ""))
        end = int(end_str.replace(",", ""))
    except (AttributeError, ValueError) as exc:
        raise ValueError(f"Invalid region format: {region}") from exc

    if not chrom or start <= 0 or end <= 0:
        raise ValueError(f"Invalid region values: {region}")
    if start > end:
        start, end = end, start
    return chrom, start, end


def to_api_region(chrom: str, start: int, end: int) -> str:
    chrom_no_prefix = chrom.replace("chr", "", 1)
    return f"{chrom_no_prefix}:{start}-{end}"


def normalize_trait_name(association: dict) -> str:
    if not isinstance(association, dict):
        return ""

    trait_name = association.get("traitName")
    if isinstance(trait_name, list):
        values = [str(value).strip() for value in trait_name if str(value).strip()]
        if values:
            return "; ".join(dict.fromkeys(values))
    elif isinstance(trait_name, str) and trait_name.strip():
        return trait_name.strip()

    efo_traits = association.get("efoTraits", [])
    labels = []
    if isinstance(efo_traits, list):
        for trait in efo_traits:
            if isinstance(trait, dict):
                label = str(trait.get("label", "")).strip()
                if label:
                    labels.append(label)
    return "; ".join(dict.fromkeys(labels))


def extract_risk_alleles(association: dict) -> list[dict]:
    if not isinstance(association, dict):
        return []

    candidates = []
    for key in ("riskAllele", "riskAlleles", "strongestRiskAlleles"):
        value = association.get(key)
        if value is None:
            continue
        if isinstance(value, list):
            candidates.extend(value)
        else:
            candidates.append(value)

    loci = association.get("loci")
    if isinstance(loci, list):
        for locus in loci:
            if not isinstance(locus, dict):
                continue
            for key in ("riskAllele", "riskAlleles", "strongestRiskAlleles"):
                value = locus.get(key)
                if value is None:
                    continue
                if isinstance(value, list):
                    candidates.extend(value)
                else:
                    candidates.append(value)

    normalized = []
    for candidate in candidates:
        if isinstance(candidate, dict):
            variant_key = str(candidate.get("key", "")).strip()
            label = str(candidate.get("label", "")).strip() or variant_key
        else:
            variant_key = str(candidate).strip()
            label = variant_key
        if variant_key.startswith("rs"):
            normalized.append({"key": variant_key, "label": label})

    unique = []
    seen = set()
    for item in normalized:
        marker = (item["key"], item["label"])
        if marker not in seen:
            seen.add(marker)
            unique.append(item)
    return unique


def read_ld_variants(path: Path) -> set[str]:
    variants: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("chr"):
                continue
            for value in line.strip().split():
                if value.startswith("rs"):
                    variants.add(value)
    return variants


def load_tag_snps_by_chromosome(path: Path) -> dict[str, list[str]]:
    chr_to_snps: dict[str, list[str]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("chr"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            rsid = parts[3].strip()
            if not rsid.startswith("rs"):
                continue
            chr_to_snps.setdefault(chrom, []).append(rsid)
    return chr_to_snps


def extract_associations(payload: dict) -> list[dict]:
    if not isinstance(payload, dict):
        return []

    if isinstance(payload.get("_embedded"), dict):
        associations = payload["_embedded"].get("associations", [])
        if isinstance(associations, list):
            return associations
    if isinstance(payload.get("associations"), list):
        return payload["associations"]
    if isinstance(payload.get("data"), list):
        return payload["data"]
    return []


def fetch_region_associations(session: requests.Session, region: str) -> list[dict]:
    all_associations: list[dict] = []
    page = 0
    while True:
        try:
            response = session.get(
                ASSOCIATIONS_URL.format(region=region),
                params={"page": page, "size": PAGE_SIZE},
                timeout=60,
            )
            response.raise_for_status()
        except requests.RequestException as exc:
            print(f"Skipping region page {region} page={page}: {exc}")
            break

        try:
            payload = response.json()
        except ValueError:
            print(f"Skipping region page {region} page={page}: invalid JSON response")
            break

        current = extract_associations(payload)
        if not current:
            break
        all_associations.extend(current)

        page_meta = payload.get("page", {})
        total_pages = page_meta.get("totalPages")
        if isinstance(total_pages, int):
            if page + 1 >= total_pages:
                break
        elif len(current) < PAGE_SIZE:
            break
        page += 1

    return all_associations


def resolve_ld_file() -> Path:
    for candidate in (Path("ld_r2_equal_higher_0.8.txt"), Path("ld_r2_equal_higher_0.8")):
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Could not find ld_r2_equal_higher_0.8(.txt)")


def main() -> int:
    ld_file = resolve_ld_file()
    tecac_file = Path("TECAC_GWAS_index_SNPs_OCT_2025")
    if not tecac_file.exists():
        raise FileNotFoundError(f"Missing required file: {tecac_file}")

    ld_variants = read_ld_variants(ld_file)
    tag_snps_by_chromosome = load_tag_snps_by_chromosome(tecac_file)
    ld_cache: dict[str, dict[str, str]] = {}
    seen_rows = set()
    rows_written = 0

    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as csv_handle:
        writer = csv.DictWriter(csv_handle, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        csv_handle.flush()

        with requests.Session() as session:
            for region in region_list:
                try:
                    chrom, start, end = parse_region(region)
                except ValueError as exc:
                    print(f"Skipping malformed region '{region}': {exc}")
                    continue

                expanded_start = max(1, start - WINDOW_BP)
                expanded_end = end + WINDOW_BP
                print(f"{chrom}:{expanded_start:,}-{expanded_end:,}")

                api_region = to_api_region(chrom, expanded_start, expanded_end)
                associations = fetch_region_associations(session, api_region)
                if not associations:
                    continue

                tag_snps = tag_snps_by_chromosome.get(chrom, [])
                if not tag_snps:
                    continue

                for association in associations:
                    if not isinstance(association, dict):
                        continue

                    trait_name = normalize_trait_name(association)
                    for risk_allele in extract_risk_alleles(association):
                        risk_key = risk_allele.get("key", "")
                        if risk_key not in ld_variants:
                            continue

                        for tag_snp in tag_snps:
                            if tag_snp not in ld_cache:
                                try:
                                    ld_rows, ld_warnings = fetch_ldtrait_rows(
                                        [tag_snp],
                                        session=session,
                                        sleep_seconds=1.0,
                                        save_raw_json=False,
                                    )
                                except Exception as exc:
                                    print(f"Skipping LD query for {tag_snp}: unexpected error: {exc}")
                                    ld_cache[tag_snp] = {}
                                    continue

                                if ld_warnings and not ld_rows:
                                    first_warning = ld_warnings[0].get("details", "No details")
                                    print(f"LD warning for {tag_snp}: {first_warning}")

                                r2_by_variant: dict[str, str] = {}
                                for ld_row in ld_rows:
                                    if not isinstance(ld_row, dict):
                                        continue
                                    rs_number = str(ld_row.get("rs_number", "")).strip()
                                    r2_value = str(ld_row.get("r2", "")).strip()
                                    if rs_number and r2_value:
                                        r2_by_variant[rs_number] = r2_value
                                ld_cache[tag_snp] = r2_by_variant

                            r2_match = ld_cache[tag_snp].get(risk_key)
                            if not r2_match:
                                continue

                            row = {
                                "riskAllele": risk_key,
                                "label": str(risk_allele.get("label", risk_key)),
                                "traitName": trait_name,
                                "tagSNP": tag_snp,
                                "r2": r2_match,
                            }
                            row_key = tuple(row[field] for field in OUTPUT_FIELDS)
                            if row_key in seen_rows:
                                continue

                            seen_rows.add(row_key)
                            writer.writerow(row)
                            csv_handle.flush()
                            rows_written += 1

    print(f"Wrote {rows_written} rows to {OUTPUT_CSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
