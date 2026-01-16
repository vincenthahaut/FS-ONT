# 1. Reads file
if(file.exists(yaml)){
  sampleSheet <- yaml::read_yaml(yaml)
} else {
  stop("YAML file does not exist - ABORT!")
}

# 2. Extract variables and check if fits requirements

# 2.1. Output directory
output_dir <- ifelse(str_sub(sampleSheet$OUTPUT_DIR, start= -1) != "/", paste0(sampleSheet$OUTPUT_DIR, "/"), sampleSheet$OUTPUT_DIR)
# 2.2. Sample ID
output_suffix <- sampleSheet$OUTPUT_SUFFIX
# 2.3. Input FASTQ
fastq_input <- sampleSheet$FASTQ_INPUT
if(!file.exists(fastq_input)){
  stop("FASTQ_INPUT - Path does not exists")
}
# 2.4. Barcoding primer sequence (N's as BC)
barcode_primer <- ifelse(is.null(sampleSheet$BARCODE_PRIMER), "CGACGCTCTTCCGATCTNNNNNNNNNNNNNAAGCAGTGGTATCAACGCAGAGT", sampleSheet$BARCODE_PRIMER)
# 2.5. Number of threads
n_cores <- as.numeric(sampleSheet$THREADS)
if(n_cores > parallel::detectCores()){
  stop("THREADS - Too many cores requested")
}
# 2.6. Number of chunks for data split
n_chunks <- as.numeric(sampleSheet$N_CHUNKS)
# 2.7. Barcode whitelist
whitelist <- case_when(sampleSheet$WHITELIST == "" ~ system.file("extdata", "whitelist.fa", package="FSnanoporeR"), TRUE ~ sampleSheet$WHITELIST)
# 2.9. Boolean - should the chimeric reads be validated or just split ?
validate_chimeric <- sampleSheet$VALIDATE_CHIMERIC
if(!is.logical(validate_chimeric)){
  stop("VALIDATE_CHIMERIC - Not boolean")
}
# 2.10. Silencing
silent <- sampleSheet$SILENT
if(!is.logical(silent)){
  stop("SILENT - Not boolean")
}
# 2.11. Minimap2 binaries / BED / .mni references
minimap2_binaries <- sampleSheet$MINIMAP2_Binaries
minimap2_bed <- sampleSheet$MINIMAP2_BED
minimap2_ref <- sampleSheet$MINIMAP2_REF
if(!file.exists(minimap2_binaries) | !file.exists(minimap2_bed) | !file.exists(minimap2_ref)){
  stop("MINIMAP2_Binaries/MINIMAP2_BED/MINIMAP2_REF - Wrong path for one of these files")
}
# 2.12. Maximum allowed memory
max_memory <- 0.75*(as.numeric(sampleSheet$MAX_MEMORY))
if(as.numeric(str_extract(as.character(memuse::Sys.meminfo()$freeram), pattern = "\\d+")) < max_memory){
  stop("Too much RAM requested")
}
# 2.13. Minimal number of reads after demultiplexing
min_reads <- sampleSheet$MIN_READS

# 2.14. Boolean - Should the reads be trimmed ?
read_trimming <- sampleSheet$READ_TRIMMING
if(!is.logical(read_trimming)){
  stop("READ_TRIMMING - Not boolean")
}
read_trimming <- ifelse(read_trimming == TRUE, "True", read_trimming)

# 2.15. Boolean - Should an HTML report be generated ?
generate_report <- sampleSheet$GENERATE_REPORT
if(!is.logical(generate_report)){
  stop("GENERATE_REPORT - Not boolean")
}
# 2.16. Boolean - Should the chunks be combined ?
combine_files <- sampleSheet$COMBINE_FILES
if(!is.logical(combine_files)){
  stop("COMBINE_FILES - Not boolean")
}
# 2.17. Should the intermediate files be removed ?
cleanup_after <- sampleSheet$CLEANUP_AFTER
if(!is.logical(cleanup_after)){
  stop("CLEANUP_AFTER - Not boolean")
}
# 2.18. Should the output folder be erased before starting ?
cleanup_before <- sampleSheet$CLEANUP_BEFORE
if(!is.logical(cleanup_before)){
  stop("CLEANUP_BEFORE - Not boolean")
}
# 2.19. What type of UMI (NONE, MONOMER, TRIMER)
umi_type <- sampleSheet$UMI_TYPE
if(!umi_type %in% c("NONE", "MONOMER", "TRIMER")){
  stop("UMI_TYPE - Not comprised in 'NONE', 'MONOMER', 'TRIMER' ")
}
# 2.20. Isoquant binaries / fasta ref / GTF ref
isoquant <- sampleSheet$ISOQUANT
isoquant_fa <- sampleSheet$ISOQUANT_FA
isoquant_gtf <- sampleSheet$ISOQUANT_GTF
if(!file.exists(isoquant_fa) | !file.exists(isoquant_gtf)){
  stop("ISOQUANT_FA/ISOQUANT_GTF - Wrong path for one of these files")
}
if(!is.logical(isoquant)){
  stop("ISOQUANT - Not boolean")
}

# 2.21. If set, downsample the reads before using isoquant
downsample_n_reads <- ifelse(sampleSheet$DOWNSAMPLE_N_READS != "", as.numeric(sampleSheet$DOWNSAMPLE_N_READS), sampleSheet$DOWNSAMPLE_N_READS)
# 2.22. Should the read orientation be guessed ?
strand_specific <- sampleSheet$STRAND_SPECIFIC
if(!is.logical(strand_specific)){
  stop("STRAND_SPECIFIC - Not boolean")
}
strand_specific <- ifelse(strand_specific == TRUE, "True", strand_specific)

# 2.23. Should only isoquant be run ?
# Will only work if you have already preprocssed the files before - useful in combination with downsampling
run_isoquant_only <- sampleSheet$RUN_ISOQUANT_ONLY
if(!is.logical(run_isoquant_only)){
  stop("RUN_ISOQUANT_ONLY - Not boolean")
}
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


# 3. Assign variables
assign("output_dir",  output_dir , envir = .GlobalEnv)
assign("output_suffix",  output_suffix , envir = .GlobalEnv)
assign("fastq_input",  fastq_input , envir = .GlobalEnv)
assign("barcode_primer",  barcode_primer , envir = .GlobalEnv)
assign("n_cores",  n_cores , envir = .GlobalEnv)
assign("n_chunks",  n_chunks , envir = .GlobalEnv)
assign("whitelist",  whitelist , envir = .GlobalEnv)
assign("validate_chimeric",  validate_chimeric , envir = .GlobalEnv)
assign("silent",  silent , envir = .GlobalEnv)
assign("minimap2_binaries",  minimap2_binaries , envir = .GlobalEnv)
assign("minimap2_bed",  minimap2_bed , envir = .GlobalEnv)
assign("minimap2_ref",  minimap2_ref , envir = .GlobalEnv)
assign("max_memory",  max_memory , envir = .GlobalEnv)
assign("min_reads",  min_reads , envir = .GlobalEnv)
assign("read_trimming",  read_trimming , envir = .GlobalEnv)
assign("generate_report",  generate_report , envir = .GlobalEnv)
assign("combine_files",  combine_files , envir = .GlobalEnv)
assign("cleanup_after",  cleanup_after , envir = .GlobalEnv)
assign("cleanup_before",  cleanup_before , envir = .GlobalEnv)
assign("umi_type",  umi_type , envir = .GlobalEnv)
assign("isoquant",  isoquant , envir = .GlobalEnv)
assign("isoquant_fa",  isoquant_fa , envir = .GlobalEnv)
assign("isoquant_gtf",  isoquant_gtf , envir = .GlobalEnv)
assign("downsample_n_reads",  downsample_n_reads , envir = .GlobalEnv)
assign("strand_specific",  strand_specific , envir = .GlobalEnv)
assign("run_isoquant_only",  run_isoquant_only , envir = .GlobalEnv)
assign("transcript_counts_type",  transcript_counts_type , envir = .GlobalEnv)
assign("selected_bc",  selected_bc , envir = .GlobalEnv)
assign("bc_detection_window",  bc_detection_window , envir = .GlobalEnv)
assign("blast_identity",  blast_identity , envir = .GlobalEnv)
assign("pcr_adapter",  pcr_adapter , envir = .GlobalEnv)
assign("overlap_threshold",  overlap_threshold , envir = .GlobalEnv)
assign("vsearch_identity",  vsearch_identity , envir = .GlobalEnv)

# 4. Save YAML to the output directory
file.copy(yaml, output_dir)
