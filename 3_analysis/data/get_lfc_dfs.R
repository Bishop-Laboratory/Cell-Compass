library(here)

if (!file.exists(here("3_analysis", "data", "lfc_dfs.RData"))) {

  library(EnsDb.Mmusculus.v79)
  library(DESeq2)
  library(tximport)
  library(readxl)
  library(tidyverse)
  
  
  # Generate metadata
  metadata <- read_tsv(
    here("3_analysis", "data", "run_metadata.txt"),
    show_col_types = FALSE
  ) %>%
    arrange(run_id) %>%
    mutate(state = case_when(
      startsWith(sample_name, "d") ~ gsub("(d\\d+)_.+", "\\1", sample_name),
      grepl("MEF", sample_name) ~ "MEF",
      grepl("iPSC", sample_name) ~ "iPSC",
      grepl("ESC", sample_name) ~ "ESC"
    )) %>%
    mutate(factors = case_when(
      grepl("OSKM", sample_name) ~ "OSKM",
      grepl("SKM", sample_name) ~ "SKM",
      TRUE ~ "NA"
    )) %>%
    mutate(adhesion = case_when(
      grepl("Epcam", sample_name) ~ "Epcam",
      grepl("GFP", sample_name) ~ "GFP",
      TRUE ~ "NA"
    )) %>%
    mutate(rep = case_when(
      grepl("_R1", sample_name) ~ "1",
      grepl("_R2", sample_name) ~ "2"
    )) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(state = factor(state,
      levels = c("MEF", "d2", "d4", "d6", "d9", "d12", "iPSC", "ESC")
    ))


  # Get tximport Object From Reprocessed SRAs
  tx2sym <- AnnotationDbi::select(
    EnsDb.Mmusculus.v79,
    keys=keys(EnsDb.Mmusculus.v79),
    columns = c("TXNAME", "SYMBOL")
  )

  run_ids <- list.dirs(
    here("2_reprocess", "quants"),
    recursive = FALSE,
    full.names = FALSE
  )

  # Keep only run_ids that exist in supplementary data
  run_ids <- run_ids[run_ids %in% metadata$run_id]

  paths <- here("2_reprocess", "quants", run_ids, "quant.sf")
  names(paths) <- run_ids

  txi <- tximport(paths,
    type = "salmon",
    tx2gene = tx2sym %>% select(TXNAME, GENEID),
    ignoreTxVersion = TRUE
  )


  # Run DESeq2
  dds <- DESeqDataSetFromTximport(
    txi,
    metadata %>% column_to_rownames("run_id"),
    ~state
  )

  dds <- DESeq(dds)


  # Create ENSMUSG to gene symbol mapping
  musg2sym <- tx2sym %>%
    select(GENEID, SYMBOL) %>%
    distinct() %>%
    mutate(SYMBOL = toupper(SYMBOL))


  # Generate lfc_dfs and save
  res_names <- resultsNames(dds)[2:8]
  deg_coefs <- 2:8
  names(deg_coefs) <- res_names

  lfc_dfs <- lapply(deg_coefs, function(coef) {
    lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE) %>%
      as.data.frame() %>%
      drop_na() %>%
      rownames_to_column("GENEID") %>%
      left_join(musg2sym, by = "GENEID")
  })

  save(lfc_dfs, file = here("3_analysis", "data", "lfc_dfs.RData"))
  rm(dds, metadata, musg2sym, tx2sym, txi, deg_coefs, paths, res_names, run_ids)
}
