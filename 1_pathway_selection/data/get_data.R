library(here)
library(tidyverse)

if (!file.exists(here("1_pathway_selection", "data", "counts.Rdata"))) {
  url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE101548&format=file"
  file_path <- here("1_pathway_selection", "data", "GSE101548_RAW.tar")
  dir_path <- here("1_pathway_selection", "data", "raws")
  
  if (!file.exists(file_path)) {
    download.file(url, file_path)
    untar(file_path, exdir = dir_path)
  }
  
  # Preprocess data and save as R variables
  raw_files <- list.files(dir_path)
  
  df <- read.table(gzfile(here(dir_path, raw_files[1])), sep = "\t", header = FALSE)
  colnames(df) <- c("gene_id", sub("_.*$", "", raw_files[1]))
  
  for (raw_file in raw_files[2:length(raw_files)]) {
    temp_df <- read.table(gzfile(here(dir_path, raw_file)), sep = "\t", header = FALSE)
    colnames(temp_df) <- c("gene_id", sub("_.*$", "", raw_file))
    df <- df %>% full_join(temp_df, by = "gene_id")
  }
  
  counts <- df %>%
    column_to_rownames("gene_id") %>%
    as.matrix()
  counts <- counts[startsWith(rownames(counts), "ENSG"), ]
  
  # Metadata from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA394766&o=acc_s%3Aa
  # Not sure how to get this programmatically
  metadata <- read.csv(here("data", "SraRunTable.txt")) %>%
    filter(Assay.Type == "RNA-Seq") %>%
    column_to_rownames("Sample.Name")
  
  save(counts, metadata, file = here("1_pathway_selection", "data", "counts.Rdata"))
  rm(df, temp_df, dir_path, file_path, raw_file, raw_files, url)
}

