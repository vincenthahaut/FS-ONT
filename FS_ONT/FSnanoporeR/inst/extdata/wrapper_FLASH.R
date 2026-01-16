#!/usr/bin/env Rscript

# WARNING_1: This script is not optimized for speed. Ubuntu 20.04 - 32 Cores, 200 Gb RAM - 1.5 million reads takes +/- 15 minutes.
# WARNING_2: This script is under active development. The parameters are set for FLASH-seq nanopore. Any major change may cause unwanted behaviors so please double check the results.

# Note: Plan to move it to snakemake only ongoing.

########################################################

# 0. Prerequists

# 0.1. Libraries
message("0. Load libraries and parse yaml file.")
suppressWarnings(suppressPackageStartupMessages(library(FSnanoporeR)))
suppressMessages(conflicted::conflict_prefer("filter", "dplyr"))
suppressMessages(conflicted::conflict_prefer("lag", "dplyr"))

# 0.2. Load Yaml
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop(paste0("A YAML file must be provided, see \n", system.file("extdata", "intialise_FS_ONT_example.yaml", package="FSnanoporeR")), call.=FALSE)
} else {
  yaml <- args[1]
  if(file.exists(yaml)){
    sampleSheet <- yaml::read_yaml(yaml)
  } else {
    stop("YAML file does not exist - ABORT!")
  }
}

# 0.3. Define Variables
output_suffix <- sampleSheet$OUTPUT_SUFFIX
output_dir <- ifelse(str_sub(sampleSheet$OUTPUT_DIR, start= -1) != "/", paste0(sampleSheet$OUTPUT_DIR, "/"), sampleSheet$OUTPUT_DIR)
fastq_input <- sampleSheet$FASTQ_INPUT
barcode_primer <- ifelse(is.null(sampleSheet$BARCODE_PRIMER), "CAGCACCTCGACGCTCTTCCGATCTNNNNNNNNNNNNNAAGCAGTGGTATCAACGCAGAGT", sampleSheet$BARCODE_PRIMER)
n_cores <- as.numeric(sampleSheet$THREADS)
n_chunks <- as.numeric(sampleSheet$N_CHUNKS)
whitelist <- case_when(sampleSheet$WHITELIST == "" ~ system.file("extdata", "whitelist.fa", package="FSnanoporeR"), TRUE ~ sampleSheet$WHITELIST)
vsearch_identity <- 0.7
validate_chimeric <- sampleSheet$VALIDATE_CHIMERIC
silent <- sampleSheet$SILENT
minimap2_binaries <- sampleSheet$MINIMAP2_Binaries
minimap2_bed <- sampleSheet$MINIMAP2_BED
minimap2_ref <- sampleSheet$MINIMAP2_REF
max_memory <- 0.75*(as.numeric(sampleSheet$MAX_MEMORY))
bc_detection_window <- 200
blast_identity <- 0.75
pcr_adapter <- "AAGCAGTGGTATCAACGCAGAGT"
overlap_threshold <- 0.8
selected_bc <- case_when(sampleSheet$SELECTED_BC == "" ~ paste0("BC", 1:384), TRUE ~ sampleSheet$SELECTED_BC)
min_reads <- sampleSheet$MIN_READS
read_trimming <- sampleSheet$READ_TRIMMING
generate_report <- sampleSheet$GENERATE_REPORT
combine_files <- sampleSheet$COMBINE_FILES
cleanup_after <- sampleSheet$CLEANUP_AFTER
cleanup_before <- sampleSheet$CLEANUP_BEFORE
umi_type <- sampleSheet$UMI_TYPE
isoquant <- sampleSheet$ISOQUANT
isoquant_fa <- sampleSheet$ISOQUANT_FA
isoquant_gtf <- sampleSheet$ISOQUANT_GTF
downsample_n_reads <- ifelse(sampleSheet$DOWNSAMPLE_N_READS != "", as.numeric(sampleSheet$DOWNSAMPLE_N_READS), sampleSheet$DOWNSAMPLE_N_READS)
strand_specific <- sampleSheet$STRAND_SPECIFIC
run_isoquant_only <- FALSE
run_isoquant_only <- sampleSheet$RUN_ISOQUANT_ONLY
transcript_counts_type <- ifelse(is.null(sampleSheet$TRANSCRIPT_COUNTS_TYPE), "unique_only", sampleSheet$TRANSCRIPT_COUNTS_TYPE)

# 0.5. Cleanup before start
# if(any(grepl(output_suffix, list.files(paste0(output_dir)))) & cleanup_before == TRUE){unlink(paste0(output_dir, "/", output_suffix, "*"), recursive = TRUE, expand = TRUE)}

# N CORES & N RAM #######################################################

if(n_cores > parallel::detectCores()){
  stop("Too many cores requested")
}

if(as.numeric(str_extract(as.character(memuse::Sys.meminfo()$freeram), pattern = "\\d+")) < max_memory){
  stop("Too much RAM requested")
}

if(run_isoquant_only != TRUE){
  # PREPARE FASTQ #######################################################
  
  message("1. Split the FASTQ into chunks for parallel processing.")
  msg0 <- prepareFASTQ(FASTQ = fastq_input,
                       OUTPUT_DIR = output_dir,
                       N_CHUNKS = n_chunks,
                       THREADS = n_cores,
                       OUTPUT_SUFFIX = output_suffix,
                       silent = silent)
  
  
  
  # SPLIT CHIMERIC #######################################################
  
  fastq_chunks <- list.files(paste0(output_dir, "/", output_suffix, "_splitFASTQ/"), pattern = "splitFASTQ", full.names = TRUE)
  
  message("2. Detect chimeric reads and split them into segments.")
  parallel::mclapply(1:length(fastq_chunks), mc.silent = silent, mc.cores = n_cores, function(x)
    
    msg2 <- split_chimeric_reads(FASTQ = fastq_chunks[x],
                                 OUTPUT_DIR = output_dir,
                                 OUTPUT_SUFFIX = paste0(output_suffix, "_", x),
                                 PCR_ADAPTER_REV = pcr_adapter,
                                 PCR_ADAPTER_FWD = pcr_adapter,
                                 CHIMERIC_DETECTION_WINDOW = bc_detection_window,
                                 IDENTITY = blast_identity,
                                 silent = silent)
    
  )
  
  
  # MAP AND VALIDATE CHIMERIC #######################################################
  
  if(validate_chimeric == TRUE){
    
    message("3. Chimeric read mapping.")
    # Minimap2 requires +/- 10Gb of memory to load the index. Asking for more than the MAX_MEMORY will kill the process.
    # Mclapply run into some issues with parallel processing + minimap2. Switched to foreach.
    minimap2_instances <- ifelse(floor(max_memory/10) < n_cores, floor(max_memory/10), n_cores)
    cl <- parallel::makeCluster(minimap2_instances)
    doParallel::registerDoParallel(cl)
    `%dopar%` <- foreach::`%dopar%`
    foreach::foreach(x=1:length(fastq_chunks), .packages = "FSnanoporeR", .export = "run_minimap2") %dopar%
      run_minimap2(
        FASTQ = c(paste0(output_dir, "/", output_suffix, "_", x, ".chimeric_reads.fq.gz"),
                  paste0(output_dir, "/", output_suffix, "_", x, ".chimeric_segments.fq.gz")),
        OUTPUT_SUFFIX = paste0(output_suffix, "_", x),
        MINIMAP2_Binaries = minimap2_binaries,
        MINIMAP2_REF = minimap2_ref,
        MINIMAP2_BED = minimap2_bed,
        ISOQUANT_MODE = FALSE,
        CHIMERIC_ON = TRUE,
        OUTPUT_DIR = output_dir,
        THREADS = 1
      )
    parallel::stopCluster(cl)
    
    message("4. Chimeric read validation.")
    parallel::mclapply(1:length(fastq_chunks), mc.silent = silent, mc.cores = n_cores, function(x)
      
      msg3 <- chimeric_mapping_validation(
        BAM_path = paste0(output_dir, output_suffix, "_", x, "_sorted.bam"),
        CHIMERIC_INFO = paste0(output_dir, output_suffix, "_", x, "_chimericSegments.info.txt"),
        OUTPUT_DIR = paste0(output_dir),
        OUTPUT_SUFFIX = paste0(output_suffix, "_", x),
        THREADS = 1,
        OVERLAP_THRESHOLD = overlap_threshold,
        silent = silent
      )
    )
    
  } else {
    message("3-4. Skip chimeric read validation.")
  }
  
  # DETECT CELL BC #######################################################
  
  message("5. Detect barcode°1 positions.")
  parallel::mclapply(1:length(fastq_chunks), mc.silent = silent, mc.cores = n_cores, function(x)
    
    msg4 <- run_vsearch(FASTQ = c(paste0(output_dir, "/", output_suffix, "_", x, ".non_chimeric_reads.fq.gz"),
                                  paste0(output_dir,  "/", output_suffix, "_", x, ".chimeric_segments.fq.gz")),
                        OUTPUT_DIR = paste0(output_dir),
                        FORWARD_SEQ = barcode_primer,
                        REVERSE_SEQ = barcode_primer,
                        SEQ_IDENTITY = vsearch_identity,
                        OUTPUT_SUFFIX = paste0(output_suffix, "_BC1_", x),
                        MIN_SEQ_LENGTH = ceiling(0.65*nchar(barcode_primer)),
                        WINDOW = bc_detection_window,
                        THREADS = 1,
                        silent = silent)
  )
  
  
  # EXTRACT AND CORRECT CELL BARCODE #######################################################
  
  message("6. Extract barcode #1 sequences and correct them.")
  parallel::mclapply(1:length(fastq_chunks), mc.silent = silent, mc.cores = n_cores, function(x)
    
    msg5 <- extractBarcodes(VSEARCH_RES_PATH = paste0(output_dir, "/vsearch_", output_suffix, "_BC1_", x, ".txt"),
                            WHITELIST_PATH = whitelist,
                            SELECTED_BC_IDs = selected_bc,
                            PRIMER_SEQ = barcode_primer,
                            THREADS = 1,
                            OUTPUT_SUFFIX = paste0(output_suffix, "_", x),
                            OUTPUT_DIR = paste0(output_dir),
                            silent = silent)
  )
  
  
  
  # TRIM - UMI - DEMULTIPLEX #######################################################
  
  message("7. Demultiplexing.")
  # I advice against running too many cores here as it will require reading the FASTQ file in memory first
  parallel::mclapply(1:length(fastq_chunks), mc.silent = silent, mc.cores = ifelse(n_cores > ceiling(max_memory/14), ceiling(max_memory/14), n_cores), function(x)
    msg6 <- filter_trim_demultiplex(BARCODES_PATH = paste0(output_dir, "/", output_suffix, "_", x, "_detected_barcodes.txt"),
                                    FASTQ = c(paste0(output_dir, "/", output_suffix, "_", x, ".non_chimeric_reads.fq.gz"),
                                              paste0(output_dir, "/", output_suffix, "_", x, ".chimeric_segments.fq.gz")),
                                    SELECTED_BC_IDs = selected_bc,
                                    UMI_TYPE = umi_type,
                                    TRIM = read_trimming,
                                    VSEARCH_BARCODE_RESULTS = paste0(output_dir, "/vsearch_", output_suffix, "_BC1_", x, ".txt"),
                                    BANLIST = ifelse(validate_chimeric == FALSE, "", paste0(output_dir, "/", output_suffix, "_", x, "_chimeric_banList.txt")),
                                    MIN_READS = floor(as.numeric(min_reads)/length(fastq_chunks)),
                                    WINDOW = bc_detection_window,
                                    OUTPUT_SUFFIX = paste0(output_suffix, "_", x),
                                    OUTPUT_DIR = paste0(output_dir, "/", output_suffix, "_", x, "_demultiplexing_pass"),
                                    THREADS = 1,
                                    STRAND_SPECIFIC = strand_specific,
                                    silent = silent)
  )
  
  
  # COMBINE FILES - MULTITHREADING #######################################################
  
  if(combine_files == TRUE){
    message("8. Combine results from parallel processing.")
    combineFiles(OUTPUT_DIR = output_dir,
                 OUTPUT_SUFFIX = output_suffix,
                 VALIDATE_CHIMERIC = validate_chimeric,
                 UMI_TYPE = umi_type,
                 THREADS = n_cores,
                 STRAND_SPECIFIC = strand_specific,
                 silent = silent)
    
  } else {
    message("Skip (8) - combining files.")
  }
  
  # GENERATE REPORT #######################################################
  
  if(generate_report == TRUE & combine_files == TRUE){
    message("9. Generate Report")
    # Really slow on big datasets
    rmd_file <- system.file("extdata", "generate_report.Rmd", package="FSnanoporeR")
    system(intern=FALSE, command=paste0("R -e 'rmarkdown::render(", "\"", rmd_file, "\"",
                                        ", params = list(yaml=", "\"", yaml, "\"", ") ",
                                        ", output_file=",
                                        "\"", output_dir, "/", output_suffix, "_FLASHseqONT_output/", output_suffix,  "_report.html", "\"", ")'"))
  } else {
    message("Skip (9) - Generate Report")
  }
  
  # CLEANUP #######################################################
  
  if(cleanup_after == TRUE){
    message("10. Clean-up")
    unlink(paste0(output_dir, "/*splitFASTQ"), recursive = TRUE)
    unlink(paste0(output_dir, "/*_*_demultiplexing_pass"), recursive = TRUE)
    unlink(paste0(output_dir, "/*_ADAPTER_rblast_db*"))
    unlink(paste0(output_dir, "/*_vsearch_read*"))
    unlink(paste0(output_dir, "/*.fq.gz"))
    unlink(paste0(output_dir, "/*report*"))
    unlink(paste0(output_dir, "/*_sorted.bam*"))
    unlink(paste0(output_dir, "/*log.txt"))
    unlink(paste0(output_dir, "/*txt"))
    
  } else {
    message("Skip (10) - Cleanup")
  }
  
}
# ISOQUANT #######################################################

if(isoquant == TRUE){
  message("11. Run Isoquant")
  
  # Isoquant tends to open too many files
  system(intern = FALSE, command = "ulimit -n 100000")
  
  run_isoquant(downsample_n_reads = downsample_n_reads,
               output_dir = output_dir,
               output_suffix = output_suffix,
               premapping = TRUE,
               n_cores = n_cores,
               isoquant_fa = isoquant_fa,
               isoquant_gtf = isoquant_gtf,
               max_memory = max_memory,
               transcript_counts_type = transcript_counts_type,
               minimap2_binaries = minimap2_binaries,
               minimap2_ref = minimap2_ref,
               minimap2_bed = minimap2_bed)
  
  isoquant_output_folder <- paste0(output_dir, "/", output_suffix, "_FLASHseqONT_output/", output_suffix, "_isoquant",
                                   ifelse(is.numeric(downsample_n_reads), paste0("_",  downsample_n_reads, "/"), "/"),
                                   output_suffix, "/")
  
  message("12. Update Isoquant read assignments")
  
  unlink(paste0(isoquant_output_folder, "/", output_suffix, "_read_assignments.updated.txt"))
  update_read_assignments(ISOQUANT_READ_ASSIGNMENTS = paste0(isoquant_output_folder, output_suffix, ".read_assignments.tsv"),
                          mode.bam = umi_type,
                          OUTPUT_SUFFIX = output_suffix,
                          OUTPUT_DIR = isoquant_output_folder,
                          silent = silent)
  
  if(!umi_type %in% c("MONOMER", "TRIMER")){
    
    message("13. No UMI - Create new count matrices")
    isoquant_count_matrix_inconsistent(UPDATED_READ_ASSIGN = paste0(isoquant_output_folder, "/", output_suffix, "_read_assignments.updated.txt"),
                                       read_bc = paste0(output_dir, output_suffix, "_FLASHseqONT_output/", output_suffix, "_detected_barcodes.txt"),
                                       OUTPUT_SUFFIX = output_suffix,
                                       OUTPUT_DIR = isoquant_output_folder,
                                       silent = silent)
    
  } else if(umi_type %in% c("MONOMER", "TRIMER")){
    unlink(paste0(isoquant_output_folder, "/umi/"), recursive = TRUE)
    
    message("13. Convert BAM Files to UMI-tools")
    
    convert_isoquant_UMItools2(ISOQUANT_BAM_FOLDER = paste0(output_dir, "/", output_suffix, "_FLASHseqONT_output/", output_suffix, "_mapped_reads",
                                                            ifelse(is.numeric(downsample_n_reads), paste0("_", downsample_n_reads, "/"), "/")),
                               mode.bam = umi_type,
                               UPDATED_READ_ASSIGN = paste0(isoquant_output_folder, "/", output_suffix, "_read_assignments.updated.txt"),
                               OUTPUT_SUFFIX = output_suffix,
                               OUTPUT_DIR = paste0(isoquant_output_folder, "/umi/"),
                               THREADS = n_cores)
    
    message("14. Run UMI-tools")
    umi.bams <- list.files(paste0(isoquant_output_folder, "/umi/"), full.names = TRUE, pattern = ".bam$")
    if(length(umi.bams) >0 ){
      unlink(paste0(isoquant_output_folder, "/umi_tools/"))
      # Run UMI_tools
      parallel::mclapply(umi.bams, mc.cores = ceiling(max_memory/10), function(x)
        run_umi_tools(BAM = x,
                      OUTPUT_SUFFIX = paste0(output_suffix, "_", str_extract(x, "BC\\d+")),
                      OUTPUT_DIR = paste0(isoquant_output_folder, "umi_tools"),
                      UMI_TOOLS_BIN = "/home/vincent.hahaut/mambaforge/envs/FSnanopore/bin/umi_tools")
      )
      
      message("15. tidy UMI-tools count tables")
      # Aggregate results
      tidy_umi_tools(UMI_TOOLS_RESULTS_FOLDER = paste0(isoquant_output_folder, "umi_tools/"),
                     OUTPUT_DIR = isoquant_output_folder,
                     OUTPUT_SUFFIX = output_suffix)
    }
  }
}

message("FINISHED!")