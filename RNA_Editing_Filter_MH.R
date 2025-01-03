#!/usr/bin/env Rscript

# Function to check if a package is installed and install it if not
check_and_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    print(paste('Installing',pkg))
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
}

# List of required packages
required_packages <- c("tidyverse", "TxDb.Hsapiens.UCSC.hg19.knownGene")
# Install missing packages
lapply(required_packages, check_and_install)

# Suppress warnings when loading libraries
library(tidyverse) # v2.0.0
library(TxDb.Hsapiens.UCSC.hg19.knownGene) #v3.2.2


# Define a function to process the data
RNA_editing_filter1 <- function(input_file) {
  cat('filter-based-removal')
  
  df <- tryCatch({
    read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("The file is empty", "\n")
    return(data.frame())  # Return an empty data frame if reading fails
  })
  
  if (nrow(df) == 0) {
    cat("No data in the file. Skipping the rest of the processing steps.\n")
    return(data.frame())  # Return an empty data frame
  }
  
  #df <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  original_df <- df
  df = df[,c(1:4,7:10)]
  # Assign column names
  colnames(df)[1:8] <- c("CHROM", "POS", "REF", "TOTAL_CNT", "MAJOR_ALLELE", "MAJOR_CNT", 
                          "MINOR_ALLELE", "MINOR_CNT")
  
  # Calculate alternative allele
  ALT = c()
  for(i in 1:nrow(df)){
    if (df$REF[i] == df$MAJOR_ALLELE[i]){
      ALT = c(ALT,df$MINOR_ALLELE[i])
    }else{
      ALT = c(ALT,df$MAJOR_ALLELE[i])
    }
  }
  df$ALT = ALT
  
  # Step 1: Apply DARNED filter
  #RNA_EDIT_PATH='resources/'
  RNA_EDIT_PATH='resources/'
  DARNED = read.table(paste0(RNA_EDIT_PATH,'DARNED_cut.txt'),sep='\t',header = T)
  colnames(DARNED)=c('CHROM','POS','REF','ALT')
  DARNED$ALT[which(DARNED$ALT=='I')]='G'
  DARNED$ALT[which(DARNED$ALT=='U')]='T'
  #darned_ids <- paste(DARNED$`#CHROM`, DARNED$POS,sep = "_")
  darned_ids <- paste(DARNED$CHROM, DARNED$POS, DARNED$REF,DARNED$ALT,sep = "_")
  df_ids <- paste(df$CHROM, df$POS,df$REF,df$ALT, sep = "_")
  # Remove rows from gtex_call if they have the same "#CHROM" and "POS" values as those in DARNED
  df_f <- df[!(df_ids %in% darned_ids), ]
  
  #Step2: Apply REDIportal filter
  REDIportal = read.table(paste0(RNA_EDIT_PATH,'REDIportal_TABLE1_hg19_cut_s.txt'),header = T)
  colnames(REDIportal)=c('CHROM','POS','REF','ALT')
  REDIportal_ids <- paste(REDIportal$CHROM, REDIportal$POS,REDIportal$REF,REDIportal$ALT, sep = "_")
  # Create a vector of combined "#CHROM" and "POS" values from gtex_call
  df_ids2 <- paste(df_f$CHROM,df_f$POS,df_f$RE,df_f$ALT, sep = "_")
  # Remove rows from gtex_call if they have the same "#CHROM" and "POS" values as those in DARNED
  df_ff <- df_f[!(df_ids2 %in% REDIportal_ids), ]
  
  
  # Step3: Strand annotation & remove potential editing sites
  # Profile reference genome
  genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
  genes_hg19_df <- as.data.frame(genes_hg19)
  
  colnames(genes_hg19_df)[1] = 'CHROM'
  
  genes_hg19_df$CHROM <- gsub("chr", "", genes_hg19_df$CHROM)
  genes_hg19_df <- genes_hg19_df %>%
    filter(CHROM %in% as.character(1:22)) %>% mutate(CHROM=as.numeric(CHROM))
  
  genes_hg19_df$strand = as.character(genes_hg19_df$strand)
  
  df_ff$strand_annotated <- NA
  
  df_ff <- df_ff[df_ff$CHROM %in% as.character(1:22), ]
  
  for (i in 1:nrow(df_ff)) {
    # Filter genes_hg19_df by chromosome
    chr_df <- genes_hg19_df %>% filter(CHROM == df_ff$CHROM[i])
    
    # Initialize an empty vector to collect matching strands
    matching_strands <- c()
    
    # Loop through chr_df to find matching positions
    for (j in 1:nrow(chr_df)) {
      if (df_ff$POS[i] >= chr_df$start[j] & df_ff$POS[i] <= chr_df$end[j]) {
        # Append the matching strand to the vector
        matching_strands <- c(matching_strands, chr_df$strand[j])
      }
    }
    
    if (length(matching_strands) > 0) {
      df_ff$strand_annotated[i] <- paste(unique(matching_strands), collapse = ",")
    }
  }
  
  # remove A>G on transcribed and T>C on untranscribed strand
  df_fff = df_ff %>% filter(!(REF=='A'&ALT=='G'&strand_annotated=='+')) # remove A>G mutations on gene transcribed strand
  df_fff = df_fff %>% filter(!(REF=='T'&ALT=='C'&strand_annotated=='-')) # remove T>C mutations on gene untranscribed strand
  
  final_df <- original_df[paste0(original_df$V1,original_df$V2) %in% paste0(df_fff$CHROM,df_fff$POS), ]
  
  #return(df_fff[,1:9])
  return(final_df)
}

RNA_editing_filter2 <- function(df) {
  cat("full-removal")
  # Read the data from the input file without headers
  df <- tryCatch({
    read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("The file is empty", "\n")
    return(data.frame())  # Return an empty data frame if reading fails
  })
  
  if (nrow(df) == 0) {
    cat("No data in the file. Skipping the rest of the processing steps.\n")
    return(data.frame())  # Return an empty data frame
  }
  
  original_df <- df
  
  #df <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  df = df[,c(1:4,7:10)]
  # Assign column names
  colnames(df)[1:8] <- c("CHROM", "POS", "REF", "TOTAL_CNT", "MAJOR_ALLELE", "MAJOR_CNT", 
                         "MINOR_ALLELE", "MINOR_CNT")
  
  # Calculate alternative allele
  ALT = c()
  for(i in 1:nrow(df)){
    if (df$REF[i] == df$MAJOR_ALLELE[i]){
      ALT = c(ALT,df$MINOR_ALLELE[i])
    }else{
      ALT = c(ALT,df$MAJOR_ALLELE[i])
    }
  }
  df$ALT = ALT
  
  df_f = df %>% filter(!(REF=='A'&ALT=='G')) # remove all A>G mutations 
  df_f = df_f %>% filter(!(REF=='T'&ALT=='C')) # remove all T>C mutations
  df_f <- df_f[df_f$CHROM %in% as.character(1:22), ]
  
  final_df <- original_df[paste0(original_df$V1,original_df$V2) %in% paste0(df_f$CHROM,df_f$POS), ]
  
  #return(df_f)  # Example processing result
  return(final_df)
}

# RNA editing filter mode 
handle_data <- function(input_file, processing_mode) {
  # Select the processing function based on the RNA-filter mode
  if (processing_mode == "filter-based-removal") {
    processed_data <- RNA_editing_filter1(input_file)
  } else if (processing_mode == "full-removal") {
    processed_data <- RNA_editing_filter2(input_file)
  } else {
    stop("Invalid processing mode. Choose either 'mode1' or 'mode2'.")
  }
  
  return(processed_data)
}

# Input & run functions
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: RNA_Editing_Filter.R <input_file> <processing_mode: full-removal/filter-based-removal> [output_file]\n")
} else {
  input_file <- as.character(args[1])  # Explicitly ensure it's treated as a string
  processing_mode <- as.character(args[2])  # Treat mode as a string too
  output_file <- ifelse(length(args) > 2, as.character(args[3]), "sSNV_post_RNA_editing_filter.tsv")
  
  # Check if the input file exists before proceeding
  if (!file.exists(input_file)) {
    stop("Error: Input file does not exist!")
  }
  
  # Process data
  processed_data <- handle_data(input_file, processing_mode)
  
  # Write the processed data to the output file
  write.table(processed_data, output_file, sep = "\t", row.names = FALSE, col.names = F, quote = FALSE)
  cat("Processed data saved to", output_file, "\n")
}
