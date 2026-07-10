# ICE-A / BEDPE Workflow For `tiles_TECAC`

This repository now supports a BEDPE-driven interaction mode compatible with public chromatin-conformation data and ICE-A style inputs.

## Input Requirements

- BEDPE file (`.bedpe` or tab-delimited text)
- First 6 columns must be:
  - `chrom1 start1 end1 chrom2 start2 end2`
- Additional columns are allowed.

## Recommended Run Modes

1. BEDPE-only (no SCREEN API; robust when API is unavailable):

```bash
python scripts/check_topological_associate_domain.py \
  --tiles-file data/inputs/tiles_TECAC \
  --assembly grch38 \
  --tissue testis \
  --bedpe data/public/encode_testis_hic_loops/testis_hic_loops_GRCh38_merged.bedpe \
  --bedpe-only
```

2. Hybrid mode (SCREEN + BEDPE together):

```bash
python scripts/check_topological_associate_domain.py \
  --tiles-file data/inputs/tiles_TECAC \
  --assembly grch38 \
  --tissue testis \
  --bedpe data/public/encode_testis_hic_loops/testis_hic_loops_GRCh38_merged.bedpe
```

## Outputs

- `outputs/tiles_TECAC/tiles_TECAC_results.json`: full structured result per query
- `outputs/tiles_TECAC/tiles_TECAC_associated_regions.tsv`: relational table (`source_region` -> `associated_region`)
- `outputs/tiles_TECAC/tiles_TECAC_links.interact`: UCSC custom track (interaction arcs)
- `outputs/tiles_TECAC/tiles_TECAC_query_report.md`: per-query report and error breakdown

## Upload To UCSC Genome Browser

1. Open UCSC Genome Browser on assembly `hg38`.
2. Go to `My Data` -> `Custom Tracks`.
3. Upload `outputs/tiles_TECAC/tiles_TECAC_links.interact`.
4. Set display mode to `full` to visualize arcs between source and partner regions.

## Focused Luciferase / H1 Browser Session

The requested luciferase-active and PAINTOR regions are stored in `data/inputs/luciferase_hesc_h1_regions.tsv`.
To regenerate the focused H1 custom tracks directly from verified ENCODE GRCh38 BEDPE loop files:

```bash
python scripts/create_luciferase_hesc_browser_track.py
```

The default source mode is `direct`. It reads released ENCODE H1 ChIA-PET loop BEDPE files under
`data/public/encode_hesc_h1_loops/`, checks SNP coordinates from `data/inputs/ld_r2_equal_higher_0.8`,
`data/inputs/TECAC_GWAS_index_SNPs_OCT_2025`, and `data/inputs/luciferase_hesc_h1_resolved_snps.bed`,
and writes source provenance into `outputs/luciferase_hesc_h1/luciferase_hesc_h1_report.md`.

Main outputs:

- `outputs/luciferase_hesc_h1/luciferase_hesc_h1_ucsc_session.txt`: combined UCSC custom-track file for hg38.
- `outputs/luciferase_hesc_h1/luciferase_hesc_h1_regions.bed`: source region track.
- `outputs/luciferase_hesc_h1/luciferase_hesc_h1_snps.bed`: SNP marker track for SNPs found in coordinate sources.
- `outputs/luciferase_hesc_h1/luciferase_hesc_h1_links.interact`: H1 interaction arcs for requested SNPs present in released ENCODE BEDPE loop files.
- `outputs/luciferase_hesc_h1/luciferase_hesc_h1_links.tsv`: per-arc audit table with source ENCODE accession and method.
- `outputs/luciferase_hesc_h1/luciferase_hesc_h1_summary.tsv` and `outputs/luciferase_hesc_h1/luciferase_hesc_h1_report.md`: per-region counts, coordinate coverage, and source provenance.
