import_yaml <- function(yaml = NULL){
  
  # 1. Reads file
  if(file.exists(yaml)){
    sampleSheet <- yaml::read_yaml(yaml)
  } else {
    stop("YAML file does not exist - ABORT!")
  }
  
  # 2. Extract variables
  
  # 2.1. Output directory
  output_dir <- ifelse(str_sub(sampleSheet$OUTPUT_DIR, start= -1) != "/", paste0(sampleSheet$OUTPUT_DIR, "/"), sampleSheet$OUTPUT_DIR)
  # 2.2. Sample ID
  output_suffix <- sampleSheet$OUTPUT_SUFFIX
  # 2.3. Input FASTQ
  fastq_input <- sampleSheet$FASTQ_INPUT
  # 2.4. Barcoding primer sequence (N's as BC)
  barcode_primer <- ifelse(is.null(sampleSheet$BARCODE_PRIMER), "CAGCACCTCGACGCTCTTCCGATCTNNNNNNNNNNNNNAAGCAGTGGTATCAACGCAGAGT", sampleSheet$BARCODE_PRIMER)
  # 2.5. Number of threads
  n_cores <- as.numeric(sampleSheet$THREADS)
  if(n_cores > parallel::detectCores()){
    stop("Too many cores requested")
  }
  # 2.6. Number of chunks for data split
  n_chunks <- as.numeric(sampleSheet$N_CHUNKS)
  # 2.7. Barcode whitelist
  whitelist <- case_when(sampleSheet$WHITELIST == "" ~ system.file("extdata", "whitelist.fa", package="FSnanoporeR"), TRUE ~ sampleSheet$WHITELIST)
  # 2.9. Boolean - should the chimeric reads be validated or just split ? 
  validate_chimeric <- sampleSheet$VALIDATE_CHIMERIC
  # 2.10. Silencing 
  silent <- sampleSheet$SILENT
  # 2.11. Minimap2 binaries / BED / .mni references
  minimap2_binaries <- sampleSheet$MINIMAP2_Binaries
  minimap2_bed <- sampleSheet$MINIMAP2_BED
  minimap2_ref <- sampleSheet$MINIMAP2_REF
  # 2.12. Maximum allowed memory
  max_memory <- 0.75*(as.numeric(sampleSheet$MAX_MEMORY))
  if(as.numeric(str_extract(as.character(memuse::Sys.meminfo()$freeram), pattern = "\\d+")) < max_memory){
    stop("Too much RAM requested")
  }
  # 2.13. Minimal number of reads after demultiplexing
  min_reads <- sampleSheet$MIN_READS
  
  # 2.14. Boolean - Should the reads be trimmed ?
  read_trimming <- sampleSheet$READ_TRIMMING
  # 2.15. Boolean - Should an HTML report be generated ? 
  generate_report <- sampleSheet$GENERATE_REPORT
  # 2.16. Boolean - Should the chunks be combined ? 
  combine_files <- sampleSheet$COMBINE_FILES
  # 2.17. Should the intermediate files be removed ? 
  cleanup_after <- sampleSheet$CLEANUP_AFTER
  # 2.18. Should the output folder be erased before starting ? 
  cleanup_before <- sampleSheet$CLEANUP_BEFORE
  # 2.19. What type of UMI (NONE, MONOMER, TRIMER)
  umi_type <- sampleSheet$UMI_TYPE
  # 2.20. Isoquant binaries / fasta ref / GTF ref
  isoquant <- sampleSheet$ISOQUANT
  isoquant_fa <- sampleSheet$ISOQUANT_FA
  isoquant_gtf <- sampleSheet$ISOQUANT_GTF
  # 2.21. If set, downsample the reads before using isoquant
  downsample_n_reads <- ifelse(sampleSheet$DOWNSAMPLE_N_READS != "", as.numeric(sampleSheet$DOWNSAMPLE_N_READS), sampleSheet$DOWNSAMPLE_N_READS)
  # 2.22. Should the read orientation be guessed ?
  strand_specific <- sampleSheet$STRAND_SPECIFIC
  # 2.23. Should only isoquant be run ?  
  # Will only work if you have already preprocssed the files before - useful in combination with downsampling
  run_isoquant_only <- sampleSheet$RUN_ISOQUANT_ONLY
  # 2.24. How should the transcript from isoquant be counted ? (mode)
  transcript_counts_type <- ifelse(is.null(sampleSheet$TRANSCRIPT_COUNTS_TYPE), "unique_only", sampleSheet$TRANSCRIPT_COUNTS_TYPE)
  # 2.25. Selected barcodes from the whitelist
  selected_bc <- case_when(sampleSheet$SELECTED_BC == "" ~ paste0("BC", 1:384), TRUE ~ sampleSheet$SELECTED_BC)
  
  # Pre-fixed parameters - should be tested if modified !!!
  
  # 2.26. Barcode detection window (bp)
  bc_detection_window <- 200
  # 2.27. BLAST search identity
  blast_identity <- 0.75
  # 2.28. dT/TSO PCR adapter sequence
  pcr_adapter <- "AAGCAGTGGTATCAACGCAGAGT"
  # 2.29. Chimeric reads validation - overlap between reads and segments
  overlap_threshold <- 0.8
  # 2.30. Barcode primer sequence search (identity)
  vsearch_identity <- 0.7
  

  
}

