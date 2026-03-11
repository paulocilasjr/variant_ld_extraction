#!/usr/bin/env python3
"""
Resolve a SNP (rsID or chr:pos) to hg38/hg19 coordinates, look up nearby or
overlapping ENCODE SCREEN cCREs, retrieve linked genes, and optionally extract
exact partner anchors from a local BEDPE chromatin-interaction file.

Examples
--------
# rsID -> hg38 position -> nearest/overlapping cCRE -> linked genes
python snp_interaction_lookup.py rs11946500 --assembly grch38 --tissue testis

# Coordinate input instead of rsID
python snp_interaction_lookup.py chr4:103439330 --assembly grch38 --tissue testis

# Same analysis, plus exact anchors from a local loop file
python snp_interaction_lookup.py rs11946500 --assembly grch38 --tissue testis \
  --bedpe testis_loops.bedpe --json

BEDPE expectations
------------------
The BEDPE file should have at least the first 6 standard columns:
chrom1 start1 end1 chrom2 start2 end2
Extra columns are preserved in the output when possible.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import json
import math
from pathlib import Path
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from dataclasses import dataclass, asdict
from typing import Any, Dict, List, Optional, Sequence, Tuple

# Ensembl REST: documented endpoints
#   https://rest.ensembl.org/documentation/info/variation_id
# SCREEN GraphQL: documented at
#   https://weng-lab.github.io/SCREEN2.0/
# Note: SCREEN site currently warns of temporary API issues; the fallback host
# below matches the official Python example host referenced in their docs.
ENSEMBL_REST = "https://rest.ensembl.org"
SCREEN_GRAPHQL_ENDPOINTS = [
    "https://screen.api.wenglab.org/graphql",
    "https://factorbook.api.wenglab.org/graphql",
]

DEFAULT_WINDOW = 10_000
DEFAULT_TIMEOUT = 30
USER_AGENT = "snp-interaction-lookup/1.0"
DEFAULT_TILES_FILE = "tiles_TECAC"
DEFAULT_TILES_TISSUE = "testis"
DEFAULT_TILES_ASSEMBLY = "grch38"
DEFAULT_RELATIONAL_OUTPUT = "tiles_TECAC_associated_regions.tsv"
DEFAULT_INTERACT_OUTPUT = "tiles_TECAC_links.interact"
DEFAULT_REPORT_OUTPUT = "tiles_TECAC_query_report.md"


class LookupErrorRuntime(RuntimeError):
    pass


@dataclass
class SNPRecord:
    query: str
    chrom: str
    pos: int  # 1-based
    assembly: str
    rsid: Optional[str] = None
    alleles: Optional[str] = None


@dataclass
class BiosampleMatch:
    name: str
    displayname: Optional[str]
    ontology: Optional[str]
    lifeStage: Optional[str]
    sampleType: Optional[str]
    query_value: Optional[str] = None
    score: float = 0.0


@dataclass
class CCRERecord:
    accession: Optional[str]
    chrom: str
    start: int
    end: int
    group: Optional[str]
    distance_to_snp: int
    overlaps_snp: bool
    ct_specific: Optional[Dict[str, Any]] = None
    raw: Optional[Dict[str, Any]] = None


@dataclass
class LinkedGeneRecord:
    accession: str
    gene: Optional[str]
    geneid: Optional[str]
    method: Optional[str]
    assay: Optional[str]
    tissue: Optional[str]
    celltype: Optional[str]
    score: Optional[float]
    displayname: Optional[str]
    source: Optional[str]
    experiment_accession: Optional[str]
    gene_coordinates: Optional[Dict[str, Any]] = None


@dataclass
class LoopAnchorRecord:
    snp_anchor: Dict[str, Any]
    partner_anchor: Dict[str, Any]
    line_number: int
    extra_fields: List[str]


@dataclass
class TileRegion:
    line_number: int
    chrom: str
    start0: int
    end0: int

    def start_1based(self) -> int:
        return self.start0 + 1

    def end_1based(self) -> int:
        return self.end0

    def midpoint_1based(self) -> int:
        return self.start0 + ((self.end0 - self.start0) // 2) + 1


@dataclass
class RegionLinkRecord:
    source_line_number: int
    source_chrom: str
    source_start0: int
    source_end0: int
    source_region: str
    query_variant: str
    association_type: str
    target_chrom: str
    target_start0: int
    target_end0: int
    associated_region: str
    target_name: str
    ccre_accession: Optional[str] = None
    method: Optional[str] = None
    tissue: Optional[str] = None
    celltype: Optional[str] = None


@dataclass
class BedpeInteraction:
    line_number: int
    chrom1: str
    start1: int
    end1: int
    chrom2: str
    start2: int
    end2: int
    name: Optional[str] = None


def eprint(*args: Any) -> None:
    print(*args, file=sys.stderr)


def normalize_assembly(assembly: str) -> str:
    a = assembly.lower().strip()
    if a in {"grch38", "hg38"}:
        return "grch38"
    if a in {"grch37", "hg19"}:
        return "grch37"
    raise ValueError(f"Unsupported assembly: {assembly}")


def http_get_json(
    url: str,
    headers: Optional[Dict[str, str]] = None,
    timeout: int = DEFAULT_TIMEOUT,
) -> Any:
    req = urllib.request.Request(
        url,
        headers={"Accept": "application/json", "User-Agent": USER_AGENT, **(headers or {})},
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise LookupErrorRuntime(f"GET {url} failed: HTTP {exc.code}: {body[:300]}") from exc
    except urllib.error.URLError as exc:
        raise LookupErrorRuntime(f"GET {url} failed: {exc}") from exc


def http_post_json(
    url: str,
    payload: Dict[str, Any],
    headers: Optional[Dict[str, str]] = None,
    timeout: int = DEFAULT_TIMEOUT,
) -> Any:
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=data,
        headers={
            "Accept": "application/json",
            "Content-Type": "application/json",
            "User-Agent": USER_AGENT,
            **(headers or {}),
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise LookupErrorRuntime(f"POST {url} failed: HTTP {exc.code}: {body[:500]}") from exc
    except urllib.error.URLError as exc:
        raise LookupErrorRuntime(f"POST {url} failed: {exc}") from exc


def graphql_query(query: str, variables: Dict[str, Any]) -> Any:
    """
    Try SCREEN hosts in order. Raise the last error if all fail.
    """
    last_error: Optional[Exception] = None
    for endpoint in SCREEN_GRAPHQL_ENDPOINTS:
        try:
            result = http_post_json(endpoint, {"query": query, "variables": variables})
            if "errors" in result:
                raise LookupErrorRuntime(f"GraphQL errors from {endpoint}: {result['errors']}")
            return result["data"]
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            time.sleep(0.4)
    assert last_error is not None
    raise LookupErrorRuntime(f"All SCREEN GraphQL endpoints failed. Last error: {last_error}")


def parse_variant_input(value: str) -> Optional[Tuple[str, int]]:
    m = re.fullmatch(r"(chr[0-9XYM]+|[0-9XYM]+):([0-9,]+)", value.strip(), re.IGNORECASE)
    if not m:
        return None
    chrom = m.group(1)
    if not chrom.lower().startswith("chr"):
        chrom = f"chr{chrom}"
    pos = int(m.group(2).replace(",", ""))
    return chrom, pos


def resolve_variant(query: str, assembly: str) -> SNPRecord:
    """
    Accept either:
      - rsID: rs11946500
      - coordinate: chr4:103439330
    """
    parsed = parse_variant_input(query)
    if parsed:
        chrom, pos = parsed
        return SNPRecord(query=query, chrom=chrom, pos=pos, assembly=assembly)

    if not query.lower().startswith("rs"):
        raise ValueError("Variant must be an rsID (e.g. rs11946500) or chr:pos (e.g. chr4:103439330).")

    url = f"{ENSEMBL_REST}/variation/human/{urllib.parse.quote(query)}?content-type=application/json"
    data = http_get_json(url)

    mappings = data.get("mappings", [])
    target_assembly = "GRCh38" if assembly == "grch38" else "GRCh37"

    selected = None
    for mapping in mappings:
        if str(mapping.get("assembly_name", "")).upper() == target_assembly.upper():
            selected = mapping
            break

    if selected is None and len(mappings) == 1:
        selected = mappings[0]

    if selected is None:
        raise LookupErrorRuntime(f"No {target_assembly} mapping found for {query}.")

    chrom = selected.get("seq_region_name") or selected.get("chr")
    start = selected.get("start")
    allele_string = selected.get("allele_string") or data.get("ancestral_allele")
    if chrom is None or start is None:
        raise LookupErrorRuntime(f"Could not parse Ensembl mapping for {query}: {selected}")

    if not str(chrom).startswith("chr"):
        chrom = f"chr{chrom}"

    return SNPRecord(
        query=query,
        chrom=str(chrom),
        pos=int(start),
        assembly=assembly,
        rsid=data.get("name", query),
        alleles=allele_string,
    )


def get_biosamples(assembly: str) -> List[Dict[str, Any]]:
    query = """
    query Biosamples($assembly: String!) {
      ccREBiosampleQuery(assembly: $assembly) {
        biosamples {
          name
          ontology
          lifeStage
          sampleType
          displayname
        }
      }
    }
    """
    data = graphql_query(query, {"assembly": assembly})
    if not isinstance(data, dict):
        return []
    biosample_block = data.get("ccREBiosampleQuery")
    if not isinstance(biosample_block, dict):
        return []
    biosamples = biosample_block.get("biosamples")
    if isinstance(biosamples, list):
        return biosamples
    if isinstance(biosamples, dict):
        return [biosamples]
    return []


def biosample_score(tissue_query: str, sample: Dict[str, Any]) -> float:
    q = tissue_query.lower().strip()
    fields = [sample.get("displayname"), sample.get("name"), sample.get("ontology"), sample.get("sampleType")]
    hay = " | ".join(str(x) for x in fields if x)
    h = hay.lower()
    if not h:
        return 0.0
    if q == h:
        return 1.0
    if q in h:
        return 0.95

    pieces = re.split(r"[^a-z0-9]+", h)
    tokens = [p for p in pieces if p]
    if q in tokens:
        return 0.9

    best = 0.0
    for token in tokens + [h]:
        best = max(best, difflib.SequenceMatcher(None, q, token).ratio())
    return best


def match_biosamples(
    samples: Sequence[Dict[str, Any]],
    tissue_query: Optional[str],
    limit: int = 5,
) -> List[BiosampleMatch]:
    if not tissue_query:
        return []

    ranked: List[BiosampleMatch] = []
    for sample in samples:
        if not isinstance(sample, dict):
            continue
        score = biosample_score(tissue_query, sample)
        if score <= 0.2:
            continue
        ranked.append(
            BiosampleMatch(
                name=sample.get("name"),
                displayname=sample.get("displayname"),
                ontology=sample.get("ontology"),
                lifeStage=sample.get("lifeStage"),
                sampleType=sample.get("sampleType"),
                query_value=sample.get("name"),
                score=score,
            )
        )

    ranked.sort(key=lambda x: (-x.score, x.displayname or x.name))
    return ranked[:limit]


def find_ccres(
    chrom: str,
    pos: int,
    assembly: str,
    biosample_query_value: Optional[str],
    window: int,
) -> List[CCRERecord]:
    start = max(1, pos - window)
    end = pos + window

    attempts = [
        ("GenomicRangeInput", {"chromosome": chrom, "start": start, "end": end}),
        ("GenomicRangeInput", {"chrom": chrom, "start": start, "end": end}),
        ("CoordinateInput", {"chromosome": chrom, "start": start, "end": end}),
        ("CoordinateInput", {"chrom": chrom, "start": start, "end": end}),
    ]

    data: Optional[Dict[str, Any]] = None
    last_error: Optional[Exception] = None
    for coord_type, coord_payload in attempts:
        query = f"""
        query CCRESearch($assembly: String!, $coords: [{coord_type}!], $cellType: String) {{
          cCRESCREENSearch(assembly: $assembly, coordinates: $coords, cellType: $cellType) {{
            chrom
            start
            len
            pct
            enhancer_zscore
            promoter_zscore
            ctcf_zscore
            dnase_zscore
            atac_zscore
            info {{
              accession
            }}
            ctspecific {{
              ct
              ctcf_zscore
              dnase_zscore
              h3k4me3_zscore
              h3k27ac_zscore
              atac_zscore
            }}
          }}
        }}
        """
        variables = {
            "assembly": assembly,
            "coords": [coord_payload],
            "cellType": biosample_query_value,
        }
        try:
            data = graphql_query(query, variables)
            break
        except Exception as exc:  # noqa: BLE001
            last_error = exc

    if data is None:
        raise LookupErrorRuntime(f"cCRESCREENSearch failed for all coordinate schema variants. Last error: {last_error}")

    raw_rows = data.get("cCRESCREENSearch", []) if isinstance(data, dict) else []
    if isinstance(raw_rows, dict):
        rows = [raw_rows]
    elif isinstance(raw_rows, list):
        rows = raw_rows
    else:
        rows = []

    out: List[CCRERecord] = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        chrom_value = row.get("chrom")
        if not chrom_value:
            continue
        try:
            ccre_start = int(row.get("start"))
            ccre_len = int(row.get("len"))
        except (TypeError, ValueError):
            continue
        if ccre_len <= 0:
            continue
        ccre_end = ccre_start + ccre_len - 1

        overlaps = ccre_start <= pos <= ccre_end
        if overlaps:
            dist = 0
        elif pos < ccre_start:
            dist = ccre_start - pos
        else:
            dist = pos - ccre_end

        ct_specific = None
        cts = row.get("ctspecific") or []
        if not isinstance(cts, list):
            cts = []
        cts_dicts = [c for c in cts if isinstance(c, dict)]
        if biosample_query_value and cts:
            selected_ct = None
            for c in cts_dicts:
                if c.get("ct") == biosample_query_value:
                    selected_ct = c
                    break
            if selected_ct is not None:
                ct_specific = selected_ct
            elif cts_dicts:
                ct_specific = cts_dicts[0]

        out.append(
            CCRERecord(
                accession=((row.get("info") or {}).get("accession") if isinstance(row.get("info"), dict) else None),
                chrom=str(chrom_value),
                start=ccre_start,
                end=ccre_end,
                group=row.get("pct"),
                distance_to_snp=dist,
                overlaps_snp=overlaps,
                ct_specific=ct_specific,
                raw=row,
            )
        )

    out.sort(key=lambda x: (x.distance_to_snp, x.start))
    return out


def get_nearest_genes_for_ccres(accessions: Sequence[str], assembly: str) -> Dict[str, List[Dict[str, Any]]]:
    if not accessions:
        return {}

    query = """
    query NearGenes($assembly: String!, $accessions: [String!]) {
      cCRESCREENSearch(assembly: $assembly, accessions: $accessions) {
        info {
          accession
        }
        nearestgenes {
          gene
          distance
        }
      }
    }
    """
    data = graphql_query(query, {"assembly": assembly, "accessions": list(accessions)})

    mapping: Dict[str, List[Dict[str, Any]]] = {}
    raw_rows = data.get("cCRESCREENSearch", []) if isinstance(data, dict) else []
    rows = raw_rows if isinstance(raw_rows, list) else ([raw_rows] if isinstance(raw_rows, dict) else [])
    for row in rows:
        if not isinstance(row, dict):
            continue
        info = row.get("info")
        acc = (info or {}).get("accession") if isinstance(info, dict) else None
        if acc:
            genes = row.get("nearestgenes")
            mapping[acc] = genes if isinstance(genes, list) else []
    return mapping


def get_linked_genes(
    accessions: Sequence[str],
    assembly: str,
    tissue_query: Optional[str] = None,
) -> List[LinkedGeneRecord]:
    if not accessions:
        return []

    query = """
    query LinkedGenes($assembly: String!, $accessions: [String!]) {
      linkedGenesQuery(assembly: $assembly, accession: $accessions) {
        accession
        p_val
        gene
        geneid
        genetype
        method
        grnaid
        effectsize
        assay
        celltype
        experiment_accession
        tissue
        variantid
        source
        slope
        score
        displayname
      }
    }
    """
    data = graphql_query(query, {"assembly": assembly, "accessions": list(accessions)})
    raw_rows = data.get("linkedGenesQuery", []) if isinstance(data, dict) else []
    rows = raw_rows if isinstance(raw_rows, list) else ([raw_rows] if isinstance(raw_rows, dict) else [])

    def keep(row: Dict[str, Any]) -> bool:
        if not tissue_query:
            return True
        q = tissue_query.lower()
        text = " | ".join(str(row.get(k, "")) for k in ("tissue", "celltype", "displayname", "assay", "source")).lower()
        if q in text:
            return True
        score = difflib.SequenceMatcher(None, q, text).ratio()
        return score >= 0.35

    filtered = [row for row in rows if isinstance(row, dict) and keep(row)]
    if tissue_query and not filtered:
        filtered = rows

    seen = set()
    out: List[LinkedGeneRecord] = []
    for row in filtered:
        key = (
            row.get("accession"),
            row.get("gene"),
            row.get("geneid"),
            row.get("method"),
            row.get("celltype"),
            row.get("tissue"),
            row.get("experiment_accession"),
        )
        if key in seen:
            continue
        seen.add(key)
        out.append(
            LinkedGeneRecord(
                accession=row.get("accession"),
                gene=row.get("gene"),
                geneid=row.get("geneid"),
                method=row.get("method"),
                assay=row.get("assay"),
                tissue=row.get("tissue"),
                celltype=row.get("celltype"),
                score=row.get("score"),
                displayname=row.get("displayname"),
                source=row.get("source"),
                experiment_accession=row.get("experiment_accession"),
            )
        )

    def rank_key(rec: LinkedGeneRecord) -> Tuple[int, float, str]:
        score = rec.score if isinstance(rec.score, (int, float)) else float("-inf")
        has_tissue = 0 if tissue_query and (
            (rec.tissue and tissue_query.lower() in rec.tissue.lower()) or
            (rec.celltype and tissue_query.lower() in rec.celltype.lower())
        ) else 1
        return (has_tissue, -(score if not math.isnan(score) else float("-inf")), rec.gene or "")

    out.sort(key=rank_key)
    return out


def get_gene_coordinates_by_name(gene_names: Sequence[str], assembly: str) -> Dict[str, Dict[str, Any]]:
    gene_names = [g for g in gene_names if g]
    if not gene_names:
        return {}

    query = """
    query GeneCoords($assembly: String!, $names: [String!]) {
      gene(assembly: $assembly, name: $names) {
        name
        id
        coordinates {
          chromosome
          start
          end
        }
      }
    }
    """
    data = graphql_query(query, {"assembly": assembly, "names": list(sorted(set(gene_names)))})

    out: Dict[str, Dict[str, Any]] = {}
    raw_rows = data.get("gene", []) if isinstance(data, dict) else []
    rows = raw_rows if isinstance(raw_rows, list) else ([raw_rows] if isinstance(raw_rows, dict) else [])
    for row in rows:
        if not isinstance(row, dict):
            continue
        name = row.get("name")
        coords = row.get("coordinates") or {}
        if name and coords:
            out[name] = {
                "gene_id": row.get("id"),
                "chrom": coords.get("chromosome"),
                "start": coords.get("start"),
                "end": coords.get("end"),
            }
    return out


def chrom_normalize(chrom: str) -> str:
    c = chrom.strip()
    if not c.lower().startswith("chr"):
        c = f"chr{c}"
    return c.lower()


def parse_bedpe(path: str, chrom: str, pos: int) -> List[LoopAnchorRecord]:
    """
    Return BEDPE records where either anchor overlaps the SNP position.
    """
    hits: List[LoopAnchorRecord] = []
    query_chrom = chrom_normalize(chrom)

    with open(path, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 6:
                continue
            if row[1].lower() == "start1":
                continue

            try:
                c1, s1, e1, c2, s2, e2 = row[:6]
                s1 = int(float(s1))
                e1 = int(float(e1))
                s2 = int(float(s2))
                e2 = int(float(e2))
            except ValueError:
                continue

            a1_match = chrom_normalize(c1) == query_chrom and s1 <= pos <= e1
            a2_match = chrom_normalize(c2) == query_chrom and s2 <= pos <= e2
            if not (a1_match or a2_match):
                continue

            if a1_match:
                snp_anchor = {"chrom": c1, "start": s1, "end": e1}
                partner_anchor = {"chrom": c2, "start": s2, "end": e2}
            else:
                snp_anchor = {"chrom": c2, "start": s2, "end": e2}
                partner_anchor = {"chrom": c1, "start": s1, "end": e1}

            hits.append(
                LoopAnchorRecord(
                    snp_anchor=snp_anchor,
                    partner_anchor=partner_anchor,
                    line_number=line_number,
                    extra_fields=row[6:],
                )
            )

    unique = {}
    for hit in hits:
        key = (
            hit.snp_anchor["chrom"], hit.snp_anchor["start"], hit.snp_anchor["end"],
            hit.partner_anchor["chrom"], hit.partner_anchor["start"], hit.partner_anchor["end"]
        )
        unique[key] = hit

    return sorted(
        unique.values(),
        key=lambda x: (x.partner_anchor["chrom"], x.partner_anchor["start"], x.partner_anchor["end"])
    )


def intervals_overlap_0based(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    # BED/BEDPE coordinates are treated as 0-based half-open intervals.
    return start_a < end_b and start_b < end_a


def load_bedpe_interactions(path: str) -> List[BedpeInteraction]:
    interactions: List[BedpeInteraction] = []
    with open(path, "r", newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 6:
                continue
            if row[1].lower() == "start1":
                continue

            try:
                c1, s1, e1, c2, s2, e2 = row[:6]
                s1 = int(float(s1))
                e1 = int(float(e1))
                s2 = int(float(s2))
                e2 = int(float(e2))
            except ValueError:
                continue

            if e1 <= s1 or e2 <= s2:
                continue

            name = None
            if len(row) > 6:
                raw_name = str(row[6]).strip()
                if raw_name and raw_name != ".":
                    name = raw_name

            interactions.append(
                BedpeInteraction(
                    line_number=line_number,
                    chrom1=ensure_chr_prefix(c1),
                    start1=s1,
                    end1=e1,
                    chrom2=ensure_chr_prefix(c2),
                    start2=s2,
                    end2=e2,
                    name=name,
                )
            )

    if not interactions:
        raise ValueError(f"No valid BEDPE interactions found in {path}")
    return interactions


def summarize_result(result: Dict[str, Any]) -> str:
    snp = result["snp"]
    lines = []

    rsid = snp.get("rsid") or snp["query"]
    allele_text = f" ({snp['alleles']})" if snp.get("alleles") else ""
    lines.append(f"SNP: {rsid} -> {snp['chrom']}:{snp['pos']}{allele_text}")
    lines.append(f"Assembly: {snp['assembly']}")

    biosample = result.get("chosen_biosample")
    if biosample:
        display = biosample.get("displayname") or biosample.get("name")
        lines.append(f"Matched biosample: {display} [score={biosample['score']:.3f}]")
    elif result.get("biosample_candidates"):
        lines.append("No single biosample selected; showing best tissue matches only.")

    ccres = result.get("ccres", [])
    if ccres:
        lines.append("")
        lines.append("Top cCRE hits:")
        for rec in ccres[:5]:
            overlap = "overlap" if rec["overlaps_snp"] else f"distance={rec['distance_to_snp']}"
            lines.append(
                f"  - {rec.get('accession') or 'NA'} "
                f"{rec['chrom']}:{rec['start']}-{rec['end']} ({overlap})"
            )
    else:
        lines.append("")
        lines.append("No cCREs returned in the requested window.")

    linked = result.get("linked_genes", [])
    if linked:
        lines.append("")
        lines.append("Linked genes:")
        for rec in linked[:10]:
            coords = rec.get("gene_coordinates") or {}
            coord_txt = ""
            if coords.get("chrom"):
                coord_txt = f" -> {coords['chrom']}:{coords['start']}-{coords['end']}"
            meta = ", ".join(x for x in [rec.get("method"), rec.get("tissue"), rec.get("celltype")] if x)
            meta_txt = f" [{meta}]" if meta else ""
            lines.append(f"  - {rec.get('gene') or rec.get('geneid') or 'NA'}{coord_txt}{meta_txt}")
    else:
        nearest = result.get("nearest_genes", {})
        if nearest:
            lines.append("")
            lines.append("Nearest genes for candidate cCREs:")
            shown = 0
            for acc, genes in nearest.items():
                for g in genes[:3]:
                    lines.append(f"  - {acc}: {g.get('gene')} (distance={g.get('distance')})")
                    shown += 1
                    if shown >= 10:
                        break
                if shown >= 10:
                    break

    loops = result.get("exact_loop_partners", [])
    if loops:
        lines.append("")
        lines.append("Exact loop partners from BEDPE:")
        for rec in loops[:10]:
            pa = rec["partner_anchor"]
            lines.append(f"  - {pa['chrom']}:{pa['start']}-{pa['end']} (line {rec['line_number']})")
    else:
        lines.append("")
        lines.append("Exact interaction anchors were not available from SCREEN alone.")
        lines.append("Provide --bedpe with a tissue-specific loop/contact file for exact partner coordinates.")

    return "\n".join(lines)


def parse_tile_regions(path: str) -> List[TileRegion]:
    regions: List[TileRegion] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue

            parts = line.split()
            if len(parts) < 3:
                eprint(f"Skipping line {line_number}: expected at least 3 columns (chrom start end)")
                continue

            chrom = parts[0]
            if not chrom.lower().startswith("chr"):
                chrom = f"chr{chrom}"

            try:
                start0 = int(parts[1])
                end0 = int(parts[2])
            except ValueError:
                eprint(f"Skipping line {line_number}: start/end are not integers")
                continue

            if start0 < 0 or end0 <= 0:
                eprint(f"Skipping line {line_number}: invalid coordinates start={start0} end={end0}")
                continue
            if end0 < start0:
                start0, end0 = end0, start0
            if end0 == start0:
                end0 = start0 + 1

            regions.append(
                TileRegion(
                    line_number=line_number,
                    chrom=chrom,
                    start0=start0,
                    end0=end0,
                )
            )

    if not regions:
        raise ValueError(f"No valid tile rows found in {path}")
    return regions


def tile_region_to_variant(region: TileRegion, mode: str) -> str:
    if mode == "start":
        pos = region.start_1based()
    elif mode == "end":
        pos = region.end_1based()
    else:
        pos = region.midpoint_1based()
    return f"{region.chrom}:{pos}"


def ensure_chr_prefix(chrom: Any) -> str:
    c = str(chrom).strip()
    if not c.lower().startswith("chr"):
        c = f"chr{c}"
    return c


def format_region_1based(chrom: str, start0: int, end0: int) -> str:
    return f"{chrom}:{start0 + 1}-{end0}"


def region_link_key(link: RegionLinkRecord) -> Tuple[Any, ...]:
    return (
        link.source_chrom,
        link.source_start0,
        link.source_end0,
        link.target_chrom,
        link.target_start0,
        link.target_end0,
        link.target_name,
        link.association_type,
    )


def deduplicate_links(links: Sequence[RegionLinkRecord]) -> List[RegionLinkRecord]:
    unique: Dict[Tuple[Any, ...], RegionLinkRecord] = {}
    for link in links:
        unique[region_link_key(link)] = link
    return sorted(
        unique.values(),
        key=lambda x: (
            x.source_line_number,
            x.association_type,
            x.target_chrom,
            x.target_start0,
            x.target_end0,
            x.target_name,
        ),
    )


def extract_links_from_bedpe_interactions(
    tile: TileRegion,
    query_variant: str,
    interactions: Sequence[BedpeInteraction],
) -> List[RegionLinkRecord]:
    source_region = format_region_1based(tile.chrom, tile.start0, tile.end0)
    query_chrom = chrom_normalize(tile.chrom)
    links: List[RegionLinkRecord] = []

    for interaction in interactions:
        a1_match = (
            chrom_normalize(interaction.chrom1) == query_chrom
            and intervals_overlap_0based(tile.start0, tile.end0, interaction.start1, interaction.end1)
        )
        a2_match = (
            chrom_normalize(interaction.chrom2) == query_chrom
            and intervals_overlap_0based(tile.start0, tile.end0, interaction.start2, interaction.end2)
        )
        if not (a1_match or a2_match):
            continue

        partners: List[Tuple[str, int, int, str]] = []
        if a1_match:
            partners.append(
                (
                    ensure_chr_prefix(interaction.chrom2),
                    interaction.start2,
                    interaction.end2,
                    interaction.name or f"bedpe_line_{interaction.line_number}",
                )
            )
        if a2_match:
            partners.append(
                (
                    ensure_chr_prefix(interaction.chrom1),
                    interaction.start1,
                    interaction.end1,
                    interaction.name or f"bedpe_line_{interaction.line_number}",
                )
            )

        for target_chrom, target_start0, target_end0, target_name in partners:
            if target_end0 <= target_start0:
                continue
            links.append(
                RegionLinkRecord(
                    source_line_number=tile.line_number,
                    source_chrom=tile.chrom,
                    source_start0=tile.start0,
                    source_end0=tile.end0,
                    source_region=source_region,
                    query_variant=query_variant,
                    association_type="bedpe_interaction",
                    target_chrom=target_chrom,
                    target_start0=target_start0,
                    target_end0=target_end0,
                    associated_region=format_region_1based(target_chrom, target_start0, target_end0),
                    target_name=target_name,
                    method=f"BEDPE:{interaction.line_number}",
                )
            )

    return deduplicate_links(links)


def extract_links_from_analysis(
    tile: TileRegion,
    query_variant: str,
    analysis: Dict[str, Any],
) -> List[RegionLinkRecord]:
    links: List[RegionLinkRecord] = []
    seen: set[Tuple[Any, ...]] = set()

    source_region = format_region_1based(tile.chrom, tile.start0, tile.end0)

    for rec in analysis.get("linked_genes", []) or []:
        coords = rec.get("gene_coordinates") or {}
        chrom = coords.get("chrom")
        start = coords.get("start")
        end = coords.get("end")
        if chrom is None or start is None or end is None:
            continue

        try:
            target_start0 = max(0, int(start) - 1)
            target_end0 = int(end)
        except (TypeError, ValueError):
            continue
        if target_end0 <= target_start0:
            continue

        target_chrom = ensure_chr_prefix(chrom)
        target_name = str(rec.get("gene") or rec.get("geneid") or "linked_gene")
        link = RegionLinkRecord(
            source_line_number=tile.line_number,
            source_chrom=tile.chrom,
            source_start0=tile.start0,
            source_end0=tile.end0,
            source_region=source_region,
            query_variant=query_variant,
            association_type="linked_gene",
            target_chrom=target_chrom,
            target_start0=target_start0,
            target_end0=target_end0,
            associated_region=format_region_1based(target_chrom, target_start0, target_end0),
            target_name=target_name,
            ccre_accession=rec.get("accession"),
            method=rec.get("method"),
            tissue=rec.get("tissue"),
            celltype=rec.get("celltype"),
        )
        key = region_link_key(link)
        if key in seen:
            continue
        seen.add(key)
        links.append(link)

    if links:
        return links

    # Fallback when no linked gene coordinates are available: link tile to
    # candidate cCRE regions so there is still a region-to-region mapping.
    for rec in analysis.get("ccres", []) or []:
        chrom = rec.get("chrom")
        start = rec.get("start")
        end = rec.get("end")
        if chrom is None or start is None or end is None:
            continue

        try:
            ccre_start = int(start)
            ccre_end = int(end)
        except (TypeError, ValueError):
            continue
        if ccre_end <= ccre_start:
            continue

        target_chrom = ensure_chr_prefix(chrom)
        target_start0 = max(0, ccre_start - 1)
        target_end0 = ccre_end
        target_name = str(rec.get("accession") or "ccre")
        link = RegionLinkRecord(
            source_line_number=tile.line_number,
            source_chrom=tile.chrom,
            source_start0=tile.start0,
            source_end0=tile.end0,
            source_region=source_region,
            query_variant=query_variant,
            association_type="ccre_candidate",
            target_chrom=target_chrom,
            target_start0=target_start0,
            target_end0=target_end0,
            associated_region=format_region_1based(target_chrom, target_start0, target_end0),
            target_name=target_name,
            ccre_accession=rec.get("accession"),
            method=None,
            tissue=None,
            celltype=None,
        )
        key = region_link_key(link)
        if key in seen:
            continue
        seen.add(key)
        links.append(link)

    return links


def write_relational_table(path: str, links: Sequence[RegionLinkRecord]) -> None:
    header = [
        "source_line_number",
        "source_region",
        "query_variant",
        "associated_region",
        "association_type",
        "target_name",
        "ccre_accession",
        "method",
        "tissue",
        "celltype",
    ]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        for link in links:
            writer.writerow(
                [
                    link.source_line_number,
                    link.source_region,
                    link.query_variant,
                    link.associated_region,
                    link.association_type,
                    link.target_name,
                    link.ccre_accession or "",
                    link.method or "",
                    link.tissue or "",
                    link.celltype or "",
                ]
            )


def to_interact_row(link: RegionLinkRecord, index: int) -> List[str]:
    span_chrom = link.source_chrom
    span_start = link.source_start0
    span_end = link.source_end0
    if link.source_chrom == link.target_chrom:
        span_start = min(link.source_start0, link.target_start0)
        span_end = max(link.source_end0, link.target_end0)

    name = f"tile{link.source_line_number}_{link.association_type}_{index}"
    color = "0,114,178" if link.association_type == "linked_gene" else "213,94,0"
    return [
        span_chrom,
        str(span_start),
        str(span_end),
        name,
        "1000",
        "1",
        ".",
        color,
        link.source_chrom,
        str(link.source_start0),
        str(link.source_end0),
        f"tile_{link.source_line_number}",
        ".",
        link.target_chrom,
        str(link.target_start0),
        str(link.target_end0),
        link.target_name,
        ".",
    ]


def write_interact_track(
    path: str,
    links: Sequence[RegionLinkRecord],
    assembly: str,
    tissue: Optional[str],
) -> None:
    tissue_value = tissue or "unspecified"
    tissue_label = tissue_value.replace(" ", "_")
    track_header = (
        f'track type=interact name="tiles_TECAC_{tissue_label}_links" '
        f'description="tiles_TECAC associated regions ({assembly} {tissue_value})" '
        "visibility=full maxHeightPixels=256:128:64"
    )
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write(track_header + "\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for idx, link in enumerate(links, start=1):
            writer.writerow(to_interact_row(link, idx))


def _shorten(text: str, limit: int = 180) -> str:
    t = (text or "").replace("\n", " ").strip()
    if len(t) <= limit:
        return t
    return t[: limit - 3] + "..."


def write_query_report(path: str, payload: Dict[str, Any]) -> None:
    results = payload.get("results", [])
    associations = payload.get("associations", [])

    assoc_by_line: Dict[int, List[str]] = {}
    if isinstance(associations, list):
        for item in associations:
            if not isinstance(item, dict):
                continue
            line_number = item.get("source_line_number")
            if not isinstance(line_number, int):
                continue
            label = (
                f"{item.get('associated_region', 'NA')} "
                f"[{item.get('association_type', 'NA')}: {item.get('target_name', 'NA')}]"
            )
            assoc_by_line.setdefault(line_number, []).append(label)

    lines: List[str] = []
    lines.append("# tiles_TECAC Query Report")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- `tiles_file`: {payload.get('tiles_file')}")
    lines.append(f"- `assembly`: {payload.get('assembly')}")
    lines.append(f"- `tissue_query`: {payload.get('tissue_query')}")
    lines.append(f"- `window`: {payload.get('window')}")
    lines.append(f"- `bedpe`: {payload.get('bedpe')}")
    sources = payload.get("sources") if isinstance(payload.get("sources"), dict) else {}
    lines.append(f"- `screen_enabled`: {sources.get('screen_enabled')}")
    lines.append(f"- `bedpe_enabled`: {sources.get('bedpe_enabled')}")
    lines.append(f"- `tile_count`: {payload.get('tile_count')}")
    lines.append(f"- `ok_count`: {payload.get('ok_count')}")
    lines.append(f"- `partial_count`: {payload.get('partial_count', 0)}")
    lines.append(f"- `error_count`: {payload.get('error_count')}")
    lines.append(f"- `association_count`: {payload.get('association_count')}")
    lines.append("")

    error_counter: Counter[str] = Counter()
    if isinstance(results, list):
        for row in results:
            if not isinstance(row, dict):
                continue
            if row.get("status") in {"error", "partial"}:
                error_counter[str(row.get("error", "unknown error"))] += 1

    lines.append("## Error Categories")
    lines.append("")
    if not error_counter:
        lines.append("- No query errors reported.")
    else:
        for error_text, count in error_counter.most_common():
            lines.append(f"- `{count}` queries: `{_shorten(error_text, 260)}`")
    lines.append("")

    lines.append("## Per-Query Results")
    lines.append("")
    lines.append("| line | source_region | query_variant | status | assoc_count | screen_assoc | bedpe_assoc | associated_regions_preview | error |")
    lines.append("|---:|---|---|---|---:|---:|---:|---|---|")

    if isinstance(results, list):
        for row in results:
            if not isinstance(row, dict):
                continue
            line_number = row.get("line_number", "")
            tile = row.get("tile") or {}
            if isinstance(tile, dict):
                try:
                    start0 = int(tile.get("start_0based", 0))
                    end0 = int(tile.get("end_0based", 0))
                except (TypeError, ValueError):
                    start0, end0 = 0, 1
                source_region = format_region_1based(
                    str(tile.get("chrom", "NA")),
                    start0,
                    end0,
                )
            else:
                source_region = "NA"

            preview_items = assoc_by_line.get(line_number if isinstance(line_number, int) else -1, [])
            preview = "; ".join(preview_items[:3])
            if len(preview_items) > 3:
                preview = f"{preview}; ... (+{len(preview_items) - 3} more)"
            if not preview:
                preview = "none"

            status = str(row.get("status", "unknown"))
            assoc_count = row.get("association_count", 0)
            source_counts = row.get("association_sources") or {}
            screen_assoc = source_counts.get("screen", 0) if isinstance(source_counts, dict) else 0
            bedpe_assoc = source_counts.get("bedpe", 0) if isinstance(source_counts, dict) else 0
            error_text = _shorten(str(row.get("error", "")), 220)
            lines.append(
                f"| {line_number} | {source_region} | {row.get('query_variant', '')} | "
                f"{status} | {assoc_count} | {screen_assoc} | {bedpe_assoc} | {preview} | {error_text} |"
            )

    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def run_analysis(
    variant: str,
    assembly: str,
    tissue: Optional[str],
    window: int,
    bedpe: Optional[str],
    biosamples: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    assembly = normalize_assembly(assembly)
    snp = resolve_variant(variant, assembly)

    active_biosamples = biosamples if biosamples is not None else get_biosamples(assembly)
    biosample_candidates = [asdict(x) for x in match_biosamples(active_biosamples, tissue)]
    chosen_biosample = biosample_candidates[0] if biosample_candidates else None
    biosample_query_value = chosen_biosample.get("query_value") if chosen_biosample else None

    ccres = find_ccres(snp.chrom, snp.pos, assembly, biosample_query_value, window)
    accessions = [r.accession for r in ccres if r.accession][:5]

    linked_genes = get_linked_genes(accessions, assembly, tissue_query=tissue)
    gene_coords = get_gene_coordinates_by_name([x.gene for x in linked_genes if x.gene], assembly)
    for rec in linked_genes:
        if rec.gene in gene_coords:
            rec.gene_coordinates = gene_coords[rec.gene]

    nearest_genes = get_nearest_genes_for_ccres(accessions, assembly)

    exact_loop_partners: List[LoopAnchorRecord] = []
    if bedpe:
        exact_loop_partners = parse_bedpe(bedpe, snp.chrom, snp.pos)

    return {
        "snp": asdict(snp),
        "tissue_query": tissue,
        "chosen_biosample": chosen_biosample,
        "biosample_candidates": biosample_candidates,
        "ccres": [asdict(r) for r in ccres[:10]],
        "linked_genes": [asdict(r) for r in linked_genes[:20]],
        "nearest_genes": nearest_genes,
        "exact_loop_partners": [asdict(r) for r in exact_loop_partners[:50]],
        "notes": {
            "exact_partner_coordinates_require_loop_data": not bool(exact_loop_partners),
            "bedpe_used": bool(bedpe),
        },
    }


def run_tiles_analysis(
    tiles_file: str,
    assembly: str,
    tissue: Optional[str],
    window: int,
    bedpe: Optional[str],
    query_mode: str,
    use_screen: bool = True,
) -> Dict[str, Any]:
    normalized_assembly = normalize_assembly(assembly)
    regions = parse_tile_regions(tiles_file)
    biosamples = get_biosamples(normalized_assembly) if use_screen else []
    bedpe_interactions: List[BedpeInteraction] = []
    if bedpe:
        bedpe_interactions = load_bedpe_interactions(bedpe)

    rows: List[Dict[str, Any]] = []
    links: List[RegionLinkRecord] = []
    ok_rows = 0
    partial_rows = 0
    error_rows = 0

    for idx, region in enumerate(regions, start=1):
        variant = tile_region_to_variant(region, query_mode)
        eprint(f"[{idx}/{len(regions)}] Processing {variant}")
        row_payload: Dict[str, Any] = {
            "line_number": region.line_number,
            "tile": {
                "chrom": region.chrom,
                "start_0based": region.start0,
                "end_0based": region.end0,
                "start_1based": region.start_1based(),
                "end_1based": region.end_1based(),
            },
            "query_variant": variant,
            "query_mode": query_mode,
        }

        row_errors: List[str] = []
        row_links: List[RegionLinkRecord] = []
        source_counts: Dict[str, int] = {}

        if use_screen:
            try:
                row_payload["analysis"] = run_analysis(
                    variant=variant,
                    assembly=normalized_assembly,
                    tissue=tissue,
                    window=window,
                    bedpe=None,
                    biosamples=biosamples,
                )
                screen_links = extract_links_from_analysis(region, variant, row_payload["analysis"])
                row_links.extend(screen_links)
                source_counts["screen"] = len(screen_links)
            except Exception as exc:  # noqa: BLE001
                row_errors.append(f"SCREEN: {exc}")
                source_counts["screen"] = 0

        if bedpe_interactions:
            try:
                bedpe_links = extract_links_from_bedpe_interactions(region, variant, bedpe_interactions)
                row_links.extend(bedpe_links)
                source_counts["bedpe"] = len(bedpe_links)
            except Exception as exc:  # noqa: BLE001
                row_errors.append(f"BEDPE: {exc}")
                source_counts["bedpe"] = 0

        row_links = deduplicate_links(row_links)
        links.extend(row_links)
        row_payload["association_count"] = len(row_links)
        row_payload["association_sources"] = source_counts

        if row_errors and not row_links:
            error_rows += 1
            row_payload["status"] = "error"
            row_payload["error"] = " | ".join(row_errors)
        elif row_errors and row_links:
            partial_rows += 1
            row_payload["status"] = "partial"
            row_payload["error"] = " | ".join(row_errors)
        else:
            ok_rows += 1
            row_payload["status"] = "ok"

        rows.append(row_payload)

    links = deduplicate_links(links)
    return {
        "tiles_file": tiles_file,
        "tile_count": len(regions),
        "query_mode": query_mode,
        "assembly": normalized_assembly,
        "tissue_query": tissue,
        "window": window,
        "bedpe": bedpe,
        "sources": {
            "screen_enabled": use_screen,
            "bedpe_enabled": bool(bedpe_interactions),
        },
        "ok_count": ok_rows,
        "partial_count": partial_rows,
        "error_count": error_rows,
        "association_count": len(links),
        "associations": [asdict(link) for link in links],
        "results": rows,
    }


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Lookup SNP-associated cCREs, linked genes, and optional exact loop partners."
    )
    p.add_argument(
        "variant",
        nargs="?",
        help="rsID (e.g. rs11946500) or chr:pos (e.g. chr4:103439330)",
    )
    p.add_argument(
        "--assembly",
        default=DEFAULT_TILES_ASSEMBLY,
        choices=["grch38", "hg38", "grch37", "hg19"],
        help="Genome assembly (default: grch38)",
    )
    p.add_argument("--tissue", default=None, help="Tissue or cell type string to match, e.g. 'testis'")
    p.add_argument(
        "--window",
        type=int,
        default=DEFAULT_WINDOW,
        help=f"Search window around SNP in bp (default: {DEFAULT_WINDOW})",
    )
    p.add_argument(
        "--bedpe",
        default=None,
        help=(
            "Optional BEDPE chromatin-interaction file (public Hi-C/HiChIP/ChIA-PET loops). "
            "In tiles mode, used to annotate interacting partner regions."
        ),
    )
    p.add_argument(
        "--bedpe-only",
        action="store_true",
        help="Tiles mode only: skip SCREEN API and annotate using --bedpe interactions only.",
    )
    p.add_argument(
        "--tiles-file",
        default=None,
        help="Optional BED-style file with chr/start/end rows (e.g. tiles_TECAC)",
    )
    p.add_argument(
        "--tile-query-mode",
        default="midpoint",
        choices=["midpoint", "start", "end"],
        help="Coordinate used from each tile for lookup (default: midpoint)",
    )
    p.add_argument(
        "--tiles-output",
        default="tiles_TECAC_results.json",
        help="Output JSON path for --tiles-file mode. Use '-' to print to stdout.",
    )
    p.add_argument(
        "--relational-output",
        default=DEFAULT_RELATIONAL_OUTPUT,
        help="Output TSV table for source/associated region pairs in --tiles-file mode.",
    )
    p.add_argument(
        "--track-output",
        default=DEFAULT_INTERACT_OUTPUT,
        help="Output UCSC interact custom-track file in --tiles-file mode.",
    )
    p.add_argument(
        "--report-output",
        default=DEFAULT_REPORT_OUTPUT,
        help="Output Markdown report with one row per query in --tiles-file mode.",
    )
    p.add_argument("--json", action="store_true", help="Print JSON instead of a text summary")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)

    auto_tiles_mode = not args.variant and not args.tiles_file
    if auto_tiles_mode:
        args.tiles_file = DEFAULT_TILES_FILE
        if not args.tissue:
            args.tissue = DEFAULT_TILES_TISSUE
        if not Path(args.tiles_file).exists():
            eprint(f"ERROR: default tiles file not found: {args.tiles_file}")
            return 1

    if args.tiles_file:
        if args.bedpe_only and not args.bedpe:
            eprint("ERROR: --bedpe-only requires --bedpe <interactions.bedpe>.")
            return 1
        try:
            payload = run_tiles_analysis(
                tiles_file=args.tiles_file,
                assembly=args.assembly,
                tissue=args.tissue,
                window=args.window,
                bedpe=args.bedpe,
                query_mode=args.tile_query_mode,
                use_screen=not args.bedpe_only,
            )
        except Exception as exc:  # noqa: BLE001
            eprint(f"ERROR: {exc}")
            return 1

        links = [RegionLinkRecord(**item) for item in payload.get("associations", [])]
        write_relational_table(args.relational_output, links)
        write_interact_track(args.track_output, links, payload["assembly"], payload.get("tissue_query"))
        write_query_report(args.report_output, payload)

        output_json = json.dumps(payload, indent=2, sort_keys=False)
        if args.tiles_output == "-":
            print(output_json)
        else:
            output_path = Path(args.tiles_output)
            output_path.write_text(output_json, encoding="utf-8")
            print(
                f"Wrote tile analysis for {payload['tile_count']} rows "
                f"to {output_path} "
                f"(ok: {payload['ok_count']}, partial: {payload.get('partial_count', 0)}, "
                f"errors: {payload['error_count']}, associations: {payload['association_count']})."
            )
        print(f"Wrote relational table: {args.relational_output}")
        print(f"Wrote custom track (interact): {args.track_output}")
        print(f"Wrote per-query report: {args.report_output}")
        return 0

    if not args.variant:
        eprint("ERROR: variant is required unless --tiles-file is provided.")
        return 1

    try:
        result = run_analysis(
            variant=args.variant,
            assembly=args.assembly,
            tissue=args.tissue,
            window=args.window,
            bedpe=args.bedpe,
        )
    except Exception as exc:  # noqa: BLE001
        eprint(f"ERROR: {exc}")
        return 1

    if args.json:
        print(json.dumps(result, indent=2, sort_keys=False))
    else:
        print(summarize_result(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
