# Check that your working dir is correct (is there a better way to do this?)
setwd("~/Cell-Compass/data")

filename <- "GSE101548_RAW.tar"
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE101548&format=file"
dir_name <- "raws"

if (!file.exists(filename)) {
  download.file(url, filename)
  untar(filename, exdir = dir_name)
}

# Use zless in terminal to inspect files
# If the data is raw, we should see only integers with occasional zeros