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

Cluster labels are preliminary interpretations based on marker genes, not validated cell-type annotations.
