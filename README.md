# variant_ld_extraction

This repository is organized so inputs, source data, scripts, and generated outputs stay separate.

## Layout

- `scripts/`: Python entry points and analysis helpers.
- `data/inputs/`: curated project inputs, including tag SNPs, LD source SNPs, luciferase regions, resolved SNP coordinates, and provenance metadata.
- `data/public/`: downloaded public datasets such as ENCODE BEDPE loop files.
- `outputs/luciferase_hesc_h1/`: active focused H1 browser tracks, audit tables, summaries, and reports.
- `outputs/tiles_TECAC/`: broader BEDPE/SCREEN tile-interaction outputs.
- `outputs/gwas_catalog/`: GWAS Catalog scan outputs.
- `outputs/cache/hesc_h1/`: historical cached H1 association table and reports.
- `docs/`: workflow notes and usage instructions.

## Main Commands

Regenerate the focused H1 UCSC browser bundle:

```bash
python scripts/create_luciferase_hesc_browser_track.py
```

Upload this file to UCSC Genome Browser on `hg38`:

```text
outputs/luciferase_hesc_h1/luciferase_hesc_h1_ucsc_session.txt
```

Inspect per-arc provenance here:

```text
outputs/luciferase_hesc_h1/luciferase_hesc_h1_links.tsv
```

Run the broader BEDPE tile workflow:

```bash
python scripts/check_topological_associate_domain.py \
  --tiles-file data/inputs/tiles_TECAC \
  --assembly grch38 \
  --tissue testis \
  --bedpe data/public/encode_testis_hic_loops/testis_hic_loops_GRCh38_merged.bedpe \
  --bedpe-only
```

## Current Focused H1 Result

The focused H1 output is generated directly from verified ENCODE GRCh38 H1 ChIA-PET loop BEDPE files and resolved SNP coordinates. The current regenerated bundle has 20 regions, 33 SNPs with coordinates, 0 missing SNP coordinates, and 322 interaction arcs.
