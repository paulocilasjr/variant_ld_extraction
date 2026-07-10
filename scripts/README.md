# Scripts

- `create_luciferase_hesc_browser_track.py`: builds the focused H1 UCSC browser bundle from curated regions, resolved SNP coordinates, and verified ENCODE BEDPE sources.
- `check_topological_associate_domain.py`: broader SNP/tile interaction lookup using SCREEN and/or BEDPE sources.
- `query_ld_link.py`: LDlink LDtrait querying helper.
- `query_tag_snps_gwas_catalog.py`: exact tag-SNP GWAS Catalog lookup.
- `scan_tag_snp_gwas_regions.py`: regional GWAS Catalog scan around tag SNPs.
- `extract_info.py`: legacy risk-allele/tag-LD extraction workflow.
- `legacy/`: older retained script snapshots. Prefer the maintained `.py` entry points above.

Run scripts from the repository root so default `data/` and `outputs/` paths resolve correctly.
