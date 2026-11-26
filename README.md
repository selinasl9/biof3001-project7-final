# BIOF3001 Project 7 — Flaviviruses in mosquito RNA-seq data

Repository layout:

- `data/` = inputs from outside the pipeline, or manual result
- `results/<stage>/` = outputs, one folder per pipeline stage
- `scripts/<stage>/` = the scripts, one folder per pipeline stage

## Contributions

### Selina: run selection, screening, phylogeny, host identification, figures and tables

Chose and filtered the SRA runs, ran the DIAMOND screen, built the alignments and the IQ-TREE phylogenies, identified host species from COI, and produced the result tables and tree figures.

| Folder | Scripts |
| --- | --- |
| `scripts/run_selection/` | `match_runs_to_logan.py`, `filter_runs_to_rna.py` |
| `scripts/screening/` | `run_diamond_blastx.sh`, `assemble_runs_without_contigs.sh` |
| `scripts/viral_identification/` | `combine_filtered_hits.py`, `extract_viral_contigs.py` |
| `scripts/alignment/` | `run_alignment.sh`, `get_cluster_representatives.py`, `fasta_tools.py` |
| `scripts/hosts/` | `lookup_host_species.py`, `merge_host_lookup.py` |
| `scripts/coi/` | `build_coi_database.sh`, `fetch_host_coi_references.py`, `identify_host_from_coi.py`, `trim_coi_alignment.py` |
| `scripts/tables/` | `build_contig_annotation.py`, `make_result_tables.py`, `summarise_bioproject.py`, `make_china_province_table.py`, `make_prevalence_matrix.py`, `lookup_reference_country.py` |
| `scripts/figures/` | `render_virus_tree.py`, `render_virus_tree_ggtree.R`, `plot_host_vs_virus.py`, `plot_cluster_summary.py` |
| `scripts/common/` | `fasta_io.py`, `geography.py`, `mosquito.py`, `newick.py`, `palette.py`, `palette.R` |

`scripts/common/` holds code shared between stages: FASTA reading, place-name and
coordinate handling, mosquito species names, Newick parsing and distances, and the
figure colours.

### Cammy: viral identification, clustering and the diversity statistics

Ran contig download, host filtering and hit filtering across screened BioProjects, clustered the viral contigs at several identity and coverage settings, and ran the ordination, PERMANOVA, PERMDISP and Kruskal-Wallis analyses.

| Folder | Scripts |
| --- | --- |
| `scripts/viral_identification/` | `download_all_logan_contigs.sh`, `run_diamond_all.sh`, `filter_hits_three_thresholds.py`, `extract_filtered_contigs.py`, `run_host_blastn.sh`, `remove_mosquito_contigs.py`, `make_full_summary.py` |
| `scripts/viral_identification/china_supplement/` | `run_CN_PRJNA1003643_host_blastn.sh`, `clean_CN_PRJNA1003643_host_hits.py`, `extract_CN_PRJNA1003643_union.py` |
| `scripts/viral_clustering/` | `run_cdhit.sh`, `run_vsearch_70.sh`, `build_matrices.py`, `filter_matrices.py`, `annotate_viral_clusters.py` |
| `scripts/statistical_analysis/` | `common.R`, `permanova.R`, `permdisp.R`, `kruskal_dunn.R`, `plot_ordination.py`, `plot_country_pcoa.R`, `plot_country_overlap.py`, `combine_statistics.py`, `summary_tables.py` |

`scripts/statistical_analysis/` R scripts call `source('scripts/statistical_analysis/common.R')`,
they are run from the repository root.
