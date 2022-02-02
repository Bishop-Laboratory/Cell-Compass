library(here)

if (!file.exists(here("3_analysis", "data", "enrichr_res.RData"))) {

  library(tidyverse)

  
  # Load relevant data
  source(here("3_analysis", "data", "get_lfc_dfs.R"))
  load(here("3_analysis", "data", "lfc_dfs.RData"))

  
  # Set parameters to filter genelist and choose pathway database
  padj_cutoff <- 0.01
  lfc_cutoff <- 0.5
  dbs <- c("GO_Biological_Process_2021")
  
  
  # Send results to enrichr
  enrichr_res <- lapply(lfc_dfs, function(df) {
    genes <- df %>%
      filter(padj <= padj_cutoff & abs(log2FoldChange) >= lfc_cutoff) %>%
      pull(SYMBOL)
    
    return(enrichR::enrichr(genes, dbs))
  })
  
  
  # Save results
  save(enrichr_res, file = here("3_analysis", "data", "enrichr_res.RData"))
  rm(dbs, lfc_cutoff, padj_cutoff)
}
