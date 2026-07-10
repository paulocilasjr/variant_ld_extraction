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
python check_topological_associate_domain.py \
  --tiles-file tiles_TECAC \
  --assembly grch38 \
  --tissue testis \
  --bedpe /path/to/public_interactions.bedpe \
  --bedpe-only
```

2. Hybrid mode (SCREEN + BEDPE together):

```bash
python check_topological_associate_domain.py \
  --tiles-file tiles_TECAC \
  --assembly grch38 \
  --tissue testis \
  --bedpe /path/to/public_interactions.bedpe
```

## Outputs

- `tiles_TECAC_results.json`: full structured result per query
- `tiles_TECAC_associated_regions.tsv`: relational table (`source_region` -> `associated_region`)
- `tiles_TECAC_links.interact`: UCSC custom track (interaction arcs)
- `tiles_TECAC_query_report.md`: per-query report and error breakdown

## Upload To UCSC Genome Browser

1. Open UCSC Genome Browser on assembly `hg38`.
2. Go to `My Data` -> `Custom Tracks`.
3. Upload `tiles_TECAC_links.interact`.
4. Set display mode to `full` to visualize arcs between source and partner regions.

## Focused Luciferase / H1-hESC Browser Session

The requested luciferase-active and PAINTOR regions are stored in `luciferase_hesc_h1_regions.tsv`.
To regenerate the focused H1-hESC custom tracks from the cached H1-hESC association table:

```bash
python create_luciferase_hesc_browser_track.py
```

Main outputs:

- `luciferase_hesc_h1_ucsc_session.txt`: combined UCSC custom-track file for hg38.
- `luciferase_hesc_h1_regions.bed`: source region track.
- `luciferase_hesc_h1_snps.bed`: SNP marker track for SNPs found in the cached source.
- `luciferase_hesc_h1_links.interact`: H1-hESC interaction arcs for requested SNPs present in the cached association table.
- `luciferase_hesc_h1_summary.tsv` and `luciferase_hesc_h1_report.md`: per-region counts and missing cached SNPs.
