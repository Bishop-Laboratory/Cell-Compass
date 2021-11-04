################################################################################
### PACKAGES
################################################################################

# Package names
packages <- c("here", "tidyverse")

# Install standard packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}


################################################################################
### MAIN
################################################################################

url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE101548&format=file"
file_path <- here("data", "GSE101548_RAW.tar")
dir_path <- here("data", "raws")

if (!file.exists(file_path)) {
  download.file(url, file_path)
  untar(file_path, exdir = dir_path)
}

# Use zless in terminal to inspect files
# If the data is raw, we should see only integers with occasional zeros

if (!file.exists(here("data", "counts.Rdata"))) {
  raw_files <- list.files(dir_path)
  
  df <- read.table(gzfile(here(dir_path, raw_files[1])), sep = "\t", header = FALSE)
  colnames(df) <- c("gene_id", str_split(raw_files[1], "_")[[1]][1])
  
  for (raw_file in raw_files[2:length(raw_files)]) {
    temp_df <- read.table(gzfile(here(dir_path, raw_file)), sep = "\t", header = FALSE)
    colnames(temp_df) <- c("gene_id", str_split(raw_file, "_")[[1]][1])
    df <- df %>% full_join(temp_df)
  }
  
  counts <- df %>%
    column_to_rownames("gene_id") %>%
    as.matrix()
  
  save(counts, file = here("data", "counts.Rdata"))
}
