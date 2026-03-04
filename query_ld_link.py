#!/usr/bin/env python3
import argparse
import csv
import json
import random
import time
from pathlib import Path
from typing import Iterable, Sequence

import requests

# Web endpoint used by the LDtrait UI (returns JSON)
ENDPOINT = "https://ldlink.nih.gov/LDlinkRestWeb/ldtrait"

# Keep these fixed
COMMON = {
    "pop": "EUR",
    "genome_build": "grch38",
    "r2_d": "r2",
    "r2_d_threshold": "0.0",
    "window": "500000",
    "ifContinue": "Continue",
}

RESULT_FIELDS = [
    "query_variant",
    "gwas_trait",
    "pmid",
    "rs_number",
    "position",
    "alleles",
    "r2",
    "d_prime",
    "risk_allele_frequency",
    "beta_or_or",
    "effect_size_95ci",
    "p_value",
    "gwas_catalog_variant",
]

WARNING_FIELDS = ["input_snp", "variant", "position", "details"]


def _safe_idx(values: Sequence[str], idx: int) -> str:
    return values[idx] if idx < len(values) else ""


def parse_ldtrait_response(data: object, input_snp: str) -> tuple[list[dict], list[dict]]:
    rows: list[dict] = []
    warnings: list[dict] = []
    if not isinstance(data, dict):
        warnings.append(
            {
                "input_snp": input_snp,
                "variant": "",
                "position": "",
                "details": f"Unexpected LDlink payload type: {type(data).__name__}",
            }
        )
        return rows, warnings

    error_message = str(data.get("error", "")).strip()
    if error_message:
        warnings.append(
            {
                "input_snp": input_snp,
                "variant": "",
                "position": "",
                "details": error_message,
            }
        )

    if "details" not in data:
        warnings.append(
            {
                "input_snp": input_snp,
                "variant": "",
                "position": "",
                "details": "Missing 'details' field in LDlink response",
            }
        )
        return rows, warnings

    details = data.get("details", {})
    if isinstance(details, str):
        warnings.append(
            {
                "input_snp": input_snp,
                "variant": "",
                "position": "",
                "details": details,
            }
        )
        return rows, warnings

    if not isinstance(details, dict):
        warnings.append(
            {
                "input_snp": input_snp,
                "variant": "",
                "position": "",
                "details": f"Unexpected details type: {type(details).__name__}",
            }
        )
        return rows, warnings

    for query_variant, qobj in details.items():
        if query_variant == "queryWarnings":
            continue
        if isinstance(qobj, dict):
            source_rows = qobj.get("aaData", [])
        elif isinstance(qobj, list):
            source_rows = qobj
        else:
            continue

        for row in source_rows:
            if not isinstance(row, list):
                continue
            rows.append(
                {
                    "query_variant": query_variant,
                    "gwas_trait": _safe_idx(row, 0),
                    "pmid": _safe_idx(row, 1),
                    "rs_number": _safe_idx(row, 2),
                    "position": _safe_idx(row, 3),
                    "alleles": _safe_idx(row, 4),
                    "r2": _safe_idx(row, 5),
                    "d_prime": _safe_idx(row, 6),
                    "risk_allele_frequency": _safe_idx(row, 8),
                    "beta_or_or": _safe_idx(row, 9),
                    "effect_size_95ci": _safe_idx(row, 10),
                    "p_value": _safe_idx(row, 11),
                    "gwas_catalog_variant": _safe_idx(row, 12),
                }
            )

    query_warnings = details.get("queryWarnings", {})
    if isinstance(query_warnings, dict):
        warning_rows = query_warnings.get("aaData", [])
    elif isinstance(query_warnings, list):
        warning_rows = query_warnings
    else:
        warning_rows = []

    for warning in warning_rows:
        if not isinstance(warning, list):
            continue
        warnings.append(
            {
                "input_snp": input_snp,
                "variant": _safe_idx(warning, 0),
                "position": _safe_idx(warning, 1),
                "details": _safe_idx(warning, 2),
            }
        )

    return rows, warnings


def fetch_ldtrait_rows(
    variants: Iterable[str],
    *,
    session: requests.Session | None = None,
    sleep_seconds: float = 1.0,
    save_raw_json: bool = False,
    raw_json_dir: Path | None = None,
    timeout: int = 180,
) -> tuple[list[dict], list[dict]]:
    rows: list[dict] = []
    warnings: list[dict] = []
    variant_list = [variant.strip() for variant in variants if variant and variant.strip()]
    if not variant_list:
        return rows, warnings

    own_session = session is None
    active_session = session or requests.Session()

    if save_raw_json:
        if raw_json_dir is None:
            raise ValueError("raw_json_dir is required when save_raw_json is True")
        raw_json_dir.mkdir(parents=True, exist_ok=True)

    try:
        for idx, snp in enumerate(variant_list):
            payload = {
                **COMMON,
                "snps": snp,
                "reference": str(random.randint(10000, 99999)),
            }
            try:
                response = active_session.post(ENDPOINT, json=payload, timeout=timeout)
                response.raise_for_status()
            except requests.RequestException as exc:
                warnings.append(
                    {
                        "input_snp": snp,
                        "variant": "",
                        "position": "",
                        "details": f"LDlink request failed: {exc}",
                    }
                )
                continue
            try:
                data = response.json()
            except ValueError:
                warnings.append(
                    {
                        "input_snp": snp,
                        "variant": "",
                        "position": "",
                        "details": "Invalid JSON from LDlink response",
                    }
                )
                continue

            if save_raw_json and raw_json_dir is not None:
                (raw_json_dir / f"{snp}.json").write_text(
                    json.dumps(data, indent=2), encoding="utf-8"
                )

            current_rows, current_warnings = parse_ldtrait_response(data, snp)
            rows.extend(current_rows)
            warnings.extend(current_warnings)

            if sleep_seconds > 0 and idx < len(variant_list) - 1:
                # LDlink asks for sequential API usage.
                time.sleep(sleep_seconds)
    finally:
        if own_session:
            active_session.close()

    return rows, warnings


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        if rows:
            writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Query LDlink LDtrait endpoint and write flattened TSV outputs."
    )
    parser.add_argument("variants", nargs="+", help="Variant IDs (for example: rs11946500)")
    parser.add_argument(
        "--out-dir",
        default="ldtrait_out",
        help="Output directory for TSVs and optional raw JSON.",
    )
    parser.add_argument(
        "--sleep-seconds",
        type=float,
        default=1.0,
        help="Delay between sequential requests (default: 1.0).",
    )
    parser.add_argument(
        "--no-raw-json",
        action="store_true",
        help="Do not write raw JSON responses.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    raw_dir = out_dir / "raw_json"

    rows, warnings = fetch_ldtrait_rows(
        args.variants,
        sleep_seconds=args.sleep_seconds,
        save_raw_json=not args.no_raw_json,
        raw_json_dir=raw_dir,
    )
    write_tsv(out_dir / "ldtrait_results.tsv", RESULT_FIELDS, rows)
    write_tsv(out_dir / "ldtrait_warnings.tsv", WARNING_FIELDS, warnings)

    print(f"Done. Results: {out_dir / 'ldtrait_results.tsv'} ({len(rows)} rows)")
    print(f"Warnings: {out_dir / 'ldtrait_warnings.tsv'} ({len(warnings)} rows)")
    if not args.no_raw_json:
        print(f"Raw JSON: {raw_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
