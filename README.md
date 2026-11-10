# GSE213688 Analysis — Sample M1

## Requirements
- R and RStudio
- The R packages `Seurat` and `dplyr`
- GSE213688 data downloaded from GEO

## Prepare the data
Extract the M1 sample into `data/GSE213688/M1`.
The M1 folder should contain `filtered_feature_bc_matrix.h5` and a `spatial/`
folder with the spot positions, tissue image, and `scalefactors_json.json`.

## Run the analysis
1. Open `analisis_M1.R` in RStudio.
2. Update the data path in `Load10X_Spatial()` to match your M1 folder.
3. Run the full script by clicking **Source**.
4. The script saves the Seurat object and marker table to `D:/`.

## Data source
GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213688

Extract the M1 sample into data/GSE213688/M1 inside the project folder. Then run spatial_transcriptomics_analysis.R in RStudio.

## Cluster interpretation

Clusters were identified from gene-expression patterns in the M1 Visium sample.
The labels below are preliminary descriptions, not validated cell-type identities.

| Cluster | Preliminary description |
|---|---|
| 0 | ECM-rich mixed stroma |
| 1 | ECM-rich mixed stroma |
| 2 | Ig-rich mixed stroma |
| 3 | SCGB2A2-high, fibrosis-rich |
| 4 | Mixed |
| 5 | Epithelium-associated |
| 6 | Adipose-associated |

## Pathologist annotations

The GEO sample includes spot-level annotations provided by the study authors/pathologist, such as Adipose tissue, Fibrosis, and Normal Epithelium. These are tissue-region annotations and are distinct from our expression-based clusters.

We compared the two using spot-level cross-tabulations. For example, 84.8% of cluster 6 spots were annotated as Adipose tissue, while cluster 6 accounted for 12.6% of all spots annotated as Adipose tissue. This suggests cluster 6 is associated with adipose regions but does not capture all adipose spots.

Data: [GSE213688](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213688)
