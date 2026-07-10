# Luciferase H1-hESC Browser Track Report

## Summary

- `regions`: 20
- `snps_requested`: 33
- `snps_found_in_cached_source`: 30
- `snps_missing_from_cached_source`: 3
- `snps_with_cached_h1_interactions`: 20
- `h1_interactions`: 304
- `h1_association_source`: ld_r2_equal_higher_0.8_hESC_H1_associated_regions.tsv

## Outputs

- `regions_bed`: `luciferase_hesc_h1_regions.bed`
- `snps_bed`: `luciferase_hesc_h1_snps.bed`
- `links_interact`: `luciferase_hesc_h1_links.interact`
- `ucsc_session`: `luciferase_hesc_h1_ucsc_session.txt`
- `summary_tsv`: `luciferase_hesc_h1_summary.tsv`
- `report_md`: `luciferase_hesc_h1_report.md`

## Per-Region H1-hESC Interactions

| region_id | region | snps | evidence | cached_h1_interactions | snps_with_interactions | note |
|---|---|---|---|---:|---|---|
| luciferase_01 | chr3:141920302-141922312 | rs9869659 rs3835166 rs3765155 | luciferase_active | 0 | none |  |
| luciferase_02 | chr4:103346459-103348483 | rs58611096 | luciferase_active | 1 | rs58611096 |  |
| luciferase_03 | chr4:103438488-103440510 | rs1194500 rs28615487 | luciferase_active | 0 | none |  |
| luciferase_04 | chr4:103450003-103452005 | rs4699071 | luciferase_active | 0 | none |  |
| luciferase_05 | chr9:124593082-124594989 | rs35685846 | luciferase_active | 36 | rs35685846 |  |
| luciferase_06 | chr9:124620553-124622642 | rs11366819 | luciferase_active | 3 | rs11366819 |  |
| luciferase_07 | chr9:124659009-124661009 | rs6478677 | luciferase_active | 8 | rs6478677 |  |
| luciferase_08 | chr9:124695030-124697031 | rs6478680 rs573194493 | luciferase_active | 4 | rs573194493 rs6478680 |  |
| luciferase_09 | chr9:124697464-124699474 | rs4838201 | luciferase_active | 9 | rs4838201 |  |
| luciferase_10 | chr9:124733875-124735971 | rs4836988 | luciferase_active | 18 | rs4836988 |  |
| luciferase_11 | chr9:124769804-124771893 | rs2184219 | luciferase_active | 23 | rs2184219 |  |
| luciferase_12 | chr19:36330894-36332896 | rs2972654 rs2967499 rs2967500 rs2918369 rs60261567 rs2972650 | luciferase_active;kallmann_syndrome_tf | 17 | rs2967499 rs2972654 | Transcription factor related to this SNP set is involved in Kallmann Syndrome |
| luciferase_13 | chr19:36377852-36379868 | rs62112634 rs73043854 | luciferase_active | 34 | rs62112634 rs73043854 |  |
| luciferase_14 | chr19:36417767-36419761 | rs58333554 | luciferase_active | 11 | rs58333554 |  |
| luciferase_15 | chr22:41021294-41023351 | rs926914 rs71327107 rs9611460 | luciferase_active | 122 | rs71327107 rs926914 rs9611460 |  |
| luciferase_16 | chr22:41090314-41092353 | rs5995992 | luciferase_active | 14 | rs5995992 |  |
| luciferase_17 | chr22:41169600-41171531 | rs6002271 rs9611509 | luciferase_active | 4 | rs6002271 rs9611509 |  |
| paintor_01 | chr7:21504068-21506112 | rs2390539 | luciferase_active;paintor | 0 | none | PAINTOR |
| paintor_02 | chr19:36266994-36269194 | rs4806292 | luciferase_active;paintor | 0 | none | PAINTOR |
| paintor_03 | chr22:41358265-41360455 | rs5751084 | luciferase_active;paintor | 0 | none | PAINTOR |

## Missing From Cached H1 Source

- rs1194500, rs2390539, rs5751084

## UCSC Use

Upload `luciferase_hesc_h1_ucsc_session.txt` as a custom track on hg38. The combined file includes region, SNP, and H1-hESC interaction tracks.
