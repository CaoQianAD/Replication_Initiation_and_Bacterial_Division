# Batch2 proteomics reproduction for Table S4

This workflow reproduces the sample-level MaxLFQ and PPM values and the dnaA and hns `Retained after preprocessing` flags for Table S4. It also exports the full protein-group matrix for the 68 manuscript samples.

## Data availability

The Batch2 DIA-NN output file (`report.parquet`) used as input to this workflow is available from the corresponding author upon request. It contains processed identification and quantification results. The expected file's SHA-256 checksum is recorded in `inputs/source_fingerprints.json`. Bruker instrument raw data are stored separately in `.d` format.

The bundled `inputs/` directory contains the run-to-sample mapping and inclusion decisions, sample design and manuscript column names, protein-group classification and gene annotation inputs, and the ordered Table S4 protein-group selection with its existing dcw annotations.

## Run the workflow

Keep `scripts/`, `inputs/`, and `requirements.txt` together. The Parquet file may be stored anywhere.

From this directory, install the Python dependencies and run the workflow:

```bash
python -m pip install -r requirements.txt
python -B scripts/reproduce_batch2.py --report "/path/to/report.parquet" --rscript Rscript
```

Use the full path to `Rscript` if it is unavailable on `PATH`. Use `--r-library` if the R dependencies are installed in a library outside R's default search paths. By default, results are written under this directory. To write them elsewhere, pass `--output` with the desired path.

## Processing details

1. **Run selection and quantification.** Each manuscript sample contributed one retained LC-MS/MS run. HA8, the hns upshift replicate A sample at approximately 152 min, was excluded after QC review. Retained records required `Global.PG.Q.Value < 0.01` and a finite, positive `PG.MaxLFQ`. When several eligible records represented the same run and protein group, the workflow took the largest MaxLFQ value.
2. **Protein-group classification and PPM.** The default mode classifies protein groups from the Batch2 report using the bundled FASTA-derived accession mapping. For each sample, it divides protein-group MaxLFQ by the sum across all eligible target protein groups and multiplies by 1,000,000. Subsequent steps carry these PPM values forward without renormalization.
3. **Experiment-specific screening.** The initial analysis matrix requires valid quantification in all retained samples. A protein group is excluded if any replicate has both a minimum PPM < 1 and a maximum PPM < 10. Downstream analysis requires a bidirectionally unique protein-group-to-gene mapping that agrees across the MaxLFQ and PPM matrices. A protein group is then excluded if its mean PPM is < 10 in any replicate.
4. **Baseline stability and model inclusion.** Within each biological replicate, the workflow divides the pre-upshift CV by the full time-course CV. The cutoff is the greater of 1.15 and the replicate-specific 90th percentile of ratios among protein groups remaining after low-abundance filtering. A protein group is excluded if its ratio exceeds the cutoff in any replicate. For numerical stability, the mean PPM used to calculate each CV is floored at 1e-6, and 1e-8 is added to the full time-course CV before forming the ratio. The time-response models use log2 PPM and require positive PPM in all retained samples. Missing values are not imputed. A `Yes` flag indicates inclusion in the corresponding model.
5. **Table S4 export.** The workflow extracts the 76 protein groups in the bundled selection, in manuscript order. The `Is_division_genes` flag uses `Analysis_Used = TRUE` in the division-gene annotation. Existing `Is_dcw_genes` values are preserved for those 76 groups and left blank for other rows in the full table.

## Outputs

| File or directory | Contents |
| --- | --- |
| `outputs/batch2/Batch2_PPM_full.csv` | PPM values and the two preprocessing flags for 2,957 protein groups quantified in at least one of the 68 manuscript samples. |
| `outputs/batch2/Batch2_MaxLFQ_linear_full.csv` and `Batch2_MaxLFQ_log2_full.csv` | Linear and log2 MaxLFQ for the same protein groups and samples. |
| `outputs/batch2/TableS4_Protein_reproduced.csv` | The 76-protein-group manuscript subset. `NA` indicates no valid quantification. |
| `intermediate/batch2/` | PPM denominators, screening audits, and the PPM model results. |
| `outputs/validation.json` | Input checksums and output counts. The validated dnaA and hns PPM models included 1,791 and 2,190 protein groups, respectively. |


