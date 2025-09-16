#!/usr/bin/env Rscript

### updated on 2025.09.16
### By: Yuchen (Alana) Cheng

# Suppress warnings when loading libraries
library(tidyverse) # v2.0.0
library(TxDb.Hsapiens.UCSC.hg19.knownGene) # v3.2.2

# =========================
# RNA editing filter 1
# =========================
RNA_editing_filter1 <- function(input_file, action = "remove") {
  cat('Running filter-based-removal\n')
  
  df <- tryCatch({
    read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("The file is empty\n")
    return(data.frame())
  })
  
  if (nrow(df) == 0) {
    cat("No data in the file. Skipping.\n")
    return(data.frame())
  }
  
  original_df <- df
  df <- df[,c(1:4,7:10)]
  colnames(df)[1:8] <- c("CHROM", "POS", "REF", "TOTAL_CNT", 
                         "MAJOR_ALLELE", "MAJOR_CNT", 
                         "MINOR_ALLELE", "MINOR_CNT")
  
  # Determine ALT
  df$ALT <- ifelse(df$REF == df$MAJOR_ALLELE, df$MINOR_ALLELE, df$MAJOR_ALLELE)
  
  # Step 1: Apply DARNED filter
  RNA_EDIT_PATH <- 'resources/'
  DARNED <- read.table(paste0(RNA_EDIT_PATH,'DARNED_cut.txt'), sep='\t', header = TRUE)
  colnames(DARNED) <- c('CHROM','POS','REF','ALT')
  DARNED$ALT[DARNED$ALT=='I'] <- 'G'
  DARNED$ALT[DARNED$ALT=='U'] <- 'T'
  darned_ids <- paste(DARNED$CHROM, DARNED$POS, DARNED$REF, DARNED$ALT, sep = "_")
  df_ids <- paste(df$CHROM, df$POS, df$REF, df$ALT, sep = "_")
  df_f <- df[!(df_ids %in% darned_ids), ]
  
  # Step 2: Apply REDIportal filter
  REDIportal <- read.table(paste0(RNA_EDIT_PATH,'REDIportal_TABLE1_hg19_cut_s.txt'), header = TRUE)
  colnames(REDIportal) <- c('CHROM','POS','REF','ALT')
  redi_ids <- paste(REDIportal$CHROM, REDIportal$POS, REDIportal$REF, REDIportal$ALT, sep = "_")
  df_ids2 <- paste(df_f$CHROM, df_f$POS, df_f$REF, df_f$ALT, sep = "_")
  df_ff <- df_f[!(df_ids2 %in% redi_ids), ]
  
  # Step 3: Strand annotation
  genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
  genes_hg19_df <- as.data.frame(genes_hg19)
  colnames(genes_hg19_df)[1] <- 'CHROM'
  genes_hg19_df$CHROM <- gsub("chr", "", genes_hg19_df$CHROM)
  genes_hg19_df <- genes_hg19_df %>%
    filter(CHROM %in% as.character(1:22)) %>%
    mutate(CHROM=as.numeric(CHROM))
  genes_hg19_df$strand <- as.character(genes_hg19_df$strand)
  
  df_ff$strand_annotated <- NA
  df_ff <- df_ff[df_ff$CHROM %in% as.character(1:22), ]
  
  for (i in 1:nrow(df_ff)) {
    chr_df <- genes_hg19_df %>% filter(CHROM == df_ff$CHROM[i])
    matching_strands <- c()
    for (j in 1:nrow(chr_df)) {
      if (df_ff$POS[i] >= chr_df$start[j] & df_ff$POS[i] <= chr_df$end[j]) {
        matching_strands <- c(matching_strands, chr_df$strand[j])
      }
    }
    if (length(matching_strands) > 0) {
      df_ff$strand_annotated[i] <- paste(unique(matching_strands), collapse = ",")
    }
  }
  
  # Remove A>G on transcribed and T>C on untranscribed strand
  df_fff <- df_ff %>% filter(!(REF=='A'&ALT=='G'&strand_annotated=='+'))
  df_fff <- df_fff %>% filter(!(REF=='T'&ALT=='C'&strand_annotated=='-'))
  print(paste("Number of rows in df_ff:",nrow(df_ff)))
  print(paste("Number of rows in df_fff:",nrow(df_fff)))
  
  # Collect IDs
  kept_ids     <- paste(df_fff$CHROM, df_fff$POS, sep = "_")
  all_ids      <- paste(df$CHROM, df$POS, sep = "_")
  original_ids <- paste(original_df$V1, original_df$V2, sep = "_")  # assumes V1=CHROM, V2=POS

  print(paste("original_ids:", head(original_ids)))
  print(paste("kept_ids:", head(kept_ids)))
  print(paste("length original_ids:",length(original_ids)))
  print(paste("length kept_ids:",length(kept_ids)))

  # Remove or flag based on whether the position is in kept_ids
  if (action == "remove") {
    final_df <- original_df[original_ids %in% kept_ids, ]
  } else if (action == "flag") {
    original_df$RNA_EDITING_FLAG <- ifelse(original_ids %in% kept_ids, "NO", "YES")
    final_df <- original_df
  }
  
  return(final_df)
}

# =========================
# RNA editing filter 2
# =========================
RNA_editing_filter2 <- function(input_file, action = "remove") {
  cat("Running full-removal\n")
  
  df <- tryCatch({
    read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("The file is empty\n")
    return(data.frame())
  })
  
  if (nrow(df) == 0) {
    cat("No data in the file. Skipping.\n")
    return(data.frame())
  }
  
  original_df <- df
  df <- df[,c(1:4,7:10)]
  colnames(df)[1:8] <- c("CHROM", "POS", "REF", "TOTAL_CNT", 
                         "MAJOR_ALLELE", "MAJOR_CNT", 
                         "MINOR_ALLELE", "MINOR_CNT")
  
  # Compute ALT allele
  df$ALT <- ifelse(df$REF == df$MAJOR_ALLELE, df$MINOR_ALLELE, df$MAJOR_ALLELE)

  # Remove all A>G and T>C mutations
  df_f <- df %>% filter(!(REF == 'A' & ALT == 'G')) %>% filter(!(REF == 'T' & ALT == 'C'))
  df_f <- df_f[df_f$CHROM %in% as.character(1:22), ]

  # Build consistent IDs (CHROM_POS)
  kept_ids     <- paste(df_f$CHROM, df_f$POS, sep = "_")
  original_ids <- paste(original_df$V1, original_df$V2, sep = "_")  # assumes V1=CHROM, V2=POS

  print(paste("original_ids:", head(original_ids)))
  print(paste("kept_ids:", head(kept_ids)))
  print(paste("length original_ids:", length(original_ids)))
  print(paste("length kept_ids:", length(kept_ids)))

  # Remove or flag based on whether the position is in kept_ids
  if (action == "remove") {
    final_df <- original_df[original_ids %in% kept_ids, ]
  } else if (action == "flag") {
    original_df$RNA_EDITING_FLAG <- ifelse(original_ids %in% kept_ids, "NO", "YES")
    final_df <- original_df
  }

  return(final_df)
}

# =========================
# Handle data
# =========================
handle_data <- function(input_file, processing_mode, action) {
  if (processing_mode == "filter-based-removal") {
    processed_data <- RNA_editing_filter1(input_file, action)
  } else if (processing_mode == "full-removal") {
    processed_data <- RNA_editing_filter2(input_file, action)
  } else {
    stop("Invalid processing mode. Choose either 'filter-based-removal' or 'full-removal'.")
  }
  return(processed_data)
}

# =========================
# Main script
# =========================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: RNA_Editing_Filter.R <input_file> <processing_mode: full-removal/filter-based-removal> <action: remove/flag> [output_file]\n")
} else {
  input_file <- as.character(args[1])
  processing_mode <- as.character(args[2])
  action <- as.character(args[3])
  
  if (length(args) > 3) {
    output_file <- as.character(args[4])
  } else {
    input_dir <- dirname(input_file)
    output_file <- file.path(input_dir, "final.RNAedit_filtered.tsv")
  }
  
  if (!file.exists(input_file)) {
    stop("Error: Input file does not exist!")
  }
  
  processed_data <- handle_data(input_file, processing_mode, action)
  write.table(processed_data, output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat("Processed data saved to", output_file, "\n")
}
