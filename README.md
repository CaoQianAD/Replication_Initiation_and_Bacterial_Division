# Replication initiation and bacterial division: Table S4 analyses

This repository provides four scripts for comprehensive analyses of the processed DNA-seq, RNA-seq and protein abundance measurements supplied in Table S4 of the associated manuscript. The scripts reproduce key manuscript analyses and generate additional exploratory results.

The published protein subset is sufficient for the protein summary script. The unpublished full proteome and full hns transcriptome matrices are not required. Additional data can be requested from the corresponding author.

## Repository contents

Place the four scripts and Table S4 in the repository root, alongside the `Analysis_annotation_files` directory. Preserve the annotation filenames and directory name.

| File or directory | Purpose and notes|
|---|---|
| dnaAupshift_analysis.R | RNA filtering, DNA dosage correction, response identification, clustering and functional enrichment after dnaA upshift |
| Protein_mean_SD_Cal.R | PPM and relative PPM means and sample SD for the protein groups supplied in Table S4 |
| Babymachine_CCR_analysis.m | Baby-machine DNA/RNA processing, per-gene correction and fitting, CCR identification and cell-number fitting |
| Babymachine_CCR_clustering_COG.R | Clustering, heatmaps, k-means diagnostics and COG/TF enrichment of the CCR genes calculated by MATLAB |
| TableS4*.xlsx | Table S4 workbook; its filename must start with `TableS4` and end with `.xlsx` for automatic detection |
| Analysis_annotation_files| Supplied gene, COG and transcription-factor annotations |
| results| Generated analysis outputs; created automatically |
| Batch2 proteomics reproduction for Table S4| reconstructs the full protein-group abundance tables and the dnaA and hns preprocessing flags used to generate the proteomics portion of Table S4. The Batch2 DIA-NN output file (`report.parquet`) required to run the workflow is available from the corresponding author upon request. See the workflow README for processing details and instructions.|

Each script contains its own analysis functions. No additional analysis scripts are required. Run Babymachine_CCR_analysis.m in `all` or `omics` mode before running Babymachine_CCR_clustering_COG.R. The other two scripts can be run independently.

## Cluster labels in the manuscript

Cluster numbers are identifiers and do not imply a biological ranking. Clusters were relabeled for presentation in the main-text figures and Data S2. The scripts retain their original output labels:

| Manuscript/Data S2 cluster | dnaA script cluster | Baby-machine script cluster |
|---|---|---|
| 1 | 2 | 4 |
| 2 | 1 | 2 |
| 3 | 4 | 3 |
| 4 | 3 | 1 |

This relabeling does not change cluster membership. K-means numeric labels may also be permuted across software versions; interpret cluster labels together with gene membership and expression profiles.
