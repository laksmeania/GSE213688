dir.create("data/GSE213688", recursive = TRUE, showWarnings = FALSE)
 getGEOSuppFiles(
     "GSE213688",
     baseDir = "data",
     makeDirectory = TRUE
 )
untar(
  "data/GSE213688/GSE213688_RAW.tar",
  exdir = "data/GSE213688/raw"
)
files <- list.files(
  "data/GSE213688/raw",
  recursive = TRUE,
  full.names = TRUE
)
head(files, 30)
sample_dirs <- list.dirs(
  "data/GSE213688/raw",
  recursive = FALSE,
  full.names = TRUE
)
