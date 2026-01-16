#!/usr/bin/env Rscript

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
  assign("sampleSheet", yaml, envir = .GlobalEnv)
  source(system.file("extdata", "import_yaml.R", package="FSnanoporeR"))
}

# 0.3. Cleanup before start
# if(any(grepl(output_suffix, list.files(paste0(output_dir)))) & cleanup_before == TRUE){unlink(paste0(output_dir, "/", output_suffix, "*"), recursive = TRUE, expand = TRUE)}

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
  # Samtools sort is implemented and the sweet spot seems to be 4-6 threads.
  # The second best is to use a lot of cores and pass everything in one go.
  # Mclapply run into some issues with parallel processing + minimap2. Switched to foreach.
  p = proc.time()

  minimap2_instances <- ifelse(floor(max_memory/10) < floor(n_cores/5), floor(max_memory/10), floor(n_cores/5))
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
      THREADS = 5
    )
  parallel::stopCluster(cl)

  proc.time() - p
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
  df <- data.frame(read_id = character(),	TYPE = character())
  sapply(1:length(fastq_chunks), function(x) write.table(df, paste0(output_dir, "/", output_suffix, "_", x, "_chimeric_banList.txt")))
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
                          MIN_READS = floor(as.numeric(min_reads)/length(fastq_chunks)),
                          THREADS = 1,
                          OUTPUT_SUFFIX = paste0(output_suffix, "_", x),
                          OUTPUT_DIR = paste0(output_dir),
                          silent = silent)
)


# DETECT UMI / ORIENTATION #####################################################

message("7. Detect UMI / Orientation.")
if(umi_type == "MONOMER"){
  parallel::mclapply(1:length(fastq_chunks), mc.silent = silent, mc.cores = n_cores, function(x)
    get_UMI(
    FASTQ = c(paste0(output_dir, "/", output_suffix, "_", x, ".non_chimeric_reads.fq.gz"),
              paste0(output_dir, "/", output_suffix, "_", x, ".chimeric_segments.fq.gz")),
    OUTPUT_DIR = paste0(output_dir, "/", output_suffix, "_", x, "_demultiplexing_pass"),
    OUTPUT_SUFFIX = paste0(output_suffix, "_", x),
    WIN = bc_detection_window,
    silent = silent,
    dT_SEQ = "AAGCAGTGGTATCAACGCAGAGTACNNNNNNNNATACTGACGCTTTTTT",
    THREADS = 1)
  )
}

if(umi_type == "TRIMER"){
  parallel::mclapply(1:length(fastq_chunks), mc.silent = silent, mc.cores = n_cores, function(x)
    get_and_correct_UMI_trimer(
    FASTQ = c(paste0(output_dir, "/", output_suffix, "_", x, ".non_chimeric_reads.fq.gz"),
              paste0(output_dir, "/", output_suffix, "_", x, ".chimeric_segments.fq.gz")),
    OUTPUT_DIR = paste0(output_dir, "/", output_suffix, "_", x, "_demultiplexing_pass"),
    OUTPUT_SUFFIX = paste0(output_suffix, "_", x),
    TRIMERS = c("TTC", "GGT", "AAA", "CCG"),
    dT_SEQ = "AAGCAGTGGTATCAACGCAGAGTNNNNNNNNNNNNNNNNNNNNNNNNCTTTTTTTTTTTTT",
    WIN = bc_detection_window,
    THREADS = 1,
    silent = silent)
  )
}

if(strand_specific == 'True' & umi_type == "NONE"){
  parallel::mclapply(1:length(fastq_chunks), mc.silent = silent, mc.cores = n_cores, function(x)
    determine_strand_orientation(
    FASTQ = c(paste0(output_dir, "/", output_suffix, "_", x, ".non_chimeric_reads.fq.gz"),
              paste0(output_dir, "/", output_suffix, "_", x, ".chimeric_segments.fq.gz")),
    OUTPUT_DIR = paste0(output_dir, "/", output_suffix, "_", x, "_demultiplexing_pass"),
    OUTPUT_SUFFIX = paste0(output_suffix, "_", x),
    dT_SEQ = "AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTT",
    THREADS = 1,
    WINDOW = bc_detection_window,
    silent = silent)
  )
}

# TRIM - UMI - DEMULTIPLEX #####################################################

message("8. Demultiplexing.")
# I advice against running too many cores here as it will require reading the FASTQ file in memory first
python_script <- system.file("extdata", "demultiplex_fastq.py", package="FSnanoporeR")
parallel::mclapply(1:length(fastq_chunks), mc.silent = silent, mc.cores = n_cores, function(x)

  system(intern = FALSE, command = paste0(
    "python ", python_script,
    " --input_fastqs_non_chimerics " , paste0(output_dir, "/", output_suffix, "_", x, ".non_chimeric_reads.fq.gz"),
    " --input_fastqs_chimerics ", paste0(output_dir, "/", output_suffix, "_", x, ".chimeric_segments.fq.gz"),
    " --umi_vsearch ", case_when(
      umi_type %in% c("MONOMER", "TRIMER") ~ paste0(output_dir, "/", output_suffix, "_", x, "_demultiplexing_pass/vsearch_", output_suffix, "_", x, "_UMI.txt"),
      umi_type == "NONE" & strand_specific == 'True' ~ paste0(output_dir, "/", output_suffix, "_", x, "_demultiplexing_pass/vsearch_", output_suffix, "_", x, "_dT.detection.txt")),
    " --umi_detected ", case_when(
      umi_type %in% c("MONOMER", "TRIMER") ~ paste0(output_dir, "/", output_suffix, "_", x, "_demultiplexing_pass/", output_suffix, "_", x, "_detected_umi.txt"),
      umi_type == "NONE" & strand_specific == 'True' ~ paste0(output_dir, "/", output_suffix, "_", x, "_demultiplexing_pass/", output_suffix, "_", x, "_detected_dT.txt")),
    " --umi_type ", umi_type,
    " --trimming ", read_trimming,
    " --strand_orientated ", strand_specific,
    " --bc_vsearch ", paste0(output_dir, "/vsearch_", output_suffix, "_BC1_", x, ".txt"),
    " --detected_bc ", paste0(output_dir, "/", output_suffix, "_", x, "_detected_barcodes.txt"),
    " --outpath ", output_dir,
    " --outsuffix ", paste0(output_suffix, "_", x),
    " --banlist_path ", paste0(output_dir, "/", output_suffix, "_", x, "_chimeric_banList.txt")
    )
  )
)


# COMBINE FILES - MULTITHREADING #######################################################

if(combine_files == TRUE){
  message("9. Combine results from parallel processing.")
  combineFiles(OUTPUT_DIR = output_dir,
               OUTPUT_SUFFIX = output_suffix,
               VALIDATE_CHIMERIC = validate_chimeric,
               UMI_TYPE = umi_type,
               THREADS = n_cores,
               STRAND_SPECIFIC = strand_specific,
               silent = silent)

} else {
  message("Skip (9) - combining files.")
}

# GENERATE REPORT #######################################################

if(generate_report == TRUE){
    message("10. Generate Report")
    # Really slow on big datasets
    rmd_file <- system.file("extdata", "generate_report_light.Rmd", package="FSnanoporeR")
    system(intern=FALSE, command=paste0("R -e 'rmarkdown::render(", "\"", rmd_file, "\"",
                                        ", params = list(yaml=", "\"", yaml, "\"", ") ",
                                        ", output_file=",
                                        "\"", output_dir, "/", output_suffix, "_FLASHseqONT_output/", output_suffix,  "_report.html", "\"", ")'"))
} else {
  message("Skip (10) - Generate Report")
}

# CLEANUP #######################################################

if(cleanup_after == TRUE){
  message("11 Clean-up")
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
  message("Skip (11) - Cleanup")
}

}
# ISOQUANT #######################################################

if(isoquant == TRUE){
  message("12. Run Isoquant")

  # Isoquant tends to open too many files
  # NOT WORKING FROM RSCRIPT: system(intern = FALSE, command = "ulimit -n 100000")

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

  message("13. Update Isoquant read assignments")

  unlink(paste0(isoquant_output_folder, "/", output_suffix, "_read_assignments.updated.txt"))
  update_read_assignments(ISOQUANT_READ_ASSIGNMENTS = paste0(isoquant_output_folder, output_suffix, ".read_assignments.tsv"),
                          mode.bam = umi_type,
                          OUTPUT_SUFFIX = output_suffix,
                          OUTPUT_DIR = isoquant_output_folder,
                          silent = silent)

  if(!umi_type %in% c("MONOMER", "TRIMER")){

    message("14. No UMI - Create new count matrices")
    isoquant_count_matrix_inconsistent(UPDATED_READ_ASSIGN = paste0(isoquant_output_folder, "/", output_suffix, "_read_assignments.updated.txt"),
                                       read_bc = paste0(output_dir, output_suffix, "_FLASHseqONT_output/", output_suffix, "_detected_barcodes.txt"),
                                       OUTPUT_SUFFIX = output_suffix,
                                       OUTPUT_DIR = isoquant_output_folder,
                                       silent = silent)

  } else if(umi_type %in% c("MONOMER", "TRIMER")){
    unlink(paste0(isoquant_output_folder, "/umi/"), recursive = TRUE)

    message("14. Convert BAM Files to UMI-tools")

    convert_isoquant_UMItools2(ISOQUANT_BAM_FOLDER = paste0(output_dir, "/", output_suffix, "_FLASHseqONT_output/", output_suffix, "_mapped_reads",
                                                            ifelse(is.numeric(downsample_n_reads), paste0("_", downsample_n_reads, "/"), "/")),
                               mode.bam = umi_type,
                               UPDATED_READ_ASSIGN = paste0(isoquant_output_folder, "/", output_suffix, "_read_assignments.updated.txt"),
                               OUTPUT_SUFFIX = output_suffix,
                               OUTPUT_DIR = paste0(isoquant_output_folder, "/umi/"),
                               THREADS = n_cores)

    message("15. Run UMI-tools")
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

      message("16. tidy UMI-tools count tables")
      # Aggregate results
      tidy_umi_tools(UMI_TOOLS_RESULTS_FOLDER = paste0(isoquant_output_folder, "umi_tools/"),
                     OUTPUT_DIR = isoquant_output_folder,
                     OUTPUT_SUFFIX = output_suffix)
    }
  }
}

message("FINISHED!")
