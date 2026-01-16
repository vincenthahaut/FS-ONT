#' @title Combine FASTQ and Barcode Files
#'
#' @description
#' Combine results from the parallel processes before demultiplexing
#'
#' @param OUTPUT_DIR character. Output directory of the other functions.
#' @param OUTPUT_SUFFIX character. Sample ID.
#' @param VALIDATE_CHIMERIC boolean. If chimeric mapping validation was run, will return the combined files.
#' @param UMI_TYPE character. "NONE", "MONOMER" or "TRIMER".
#' @param STRAND_SPECIFIC boolean, were the strand specific files generated?
#' @param silent boolean. Should some extra info be shown ?
#' @param THREADS numeric.
#'
#' @export
combineFiles <- function(OUTPUT_DIR = NULL, OUTPUT_SUFFIX = NULL, VALIDATE_CHIMERIC = TRUE, UMI_TYPE = "NONE", THREADS = 4, silent = TRUE, STRAND_SPECIFIC = NULL){

  suppressWarnings(dir.create(paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_FLASHseqONT_output/")))

  # 1. Text files
  if(!silent){message("Text files")}

  patterns <- data.frame(
    pattern = c("_blast_search.txt",
                "_chimericSegments.info.txt",
                "_chimeric_banList.txt",
                "_detected_barcodes.txt",
                "_in_silico_chimerics.log.txt",
                "_mapping_chimerics.log.txt",
                "vsearch_*BC1*.txt",
                "vsearch_.*_UMI*.txt",
                "vsearch_*_dT.detection.txt",
                "_detected_barcodes.unfiltered.txt",
                "_extractBarcodes.report.txt",
                "_splitChimeric_reads.report.txt",
                "_chimericReads.ids.txt",
                "_demultiplex.report.txt",
                "_get_UMI.report.txt",
                "_detected_umi.txt",
                "_read_length.report.txt",
                "_detected_dT.txt",
                # FASTQ
                ".chimeric_reads.fq.gz",
                ".chimeric_segments.fq.gz",
                ".non_chimeric_reads.fq.gz",
                # BAM
                "_chimeric_sorted.bam"),
    file_suffix = c("_blast_search.txt",
                    "_chimericSegments.info.txt",
                    "_chimeric_banList.txt",
                    "_detected_barcodes.txt",
                    "_in_silico_chimerics.log.txt",
                    "_mapping_chimerics.log.txt",
                    "_vsearch.barcodes.txt",
                    "_vsearch.umi.txt",
                    "_vsearch.dT.strandness.txt",
                    "_detected_barcodes.unfiltered.txt",
                    "_extractBarcodes.report.txt",
                    "_splitChimeric_reads.report.txt",
                    "_chimericReads.ids.txt",
                    "_demultiplex.report.txt",
                    "_get_UMI.report.txt",
                    "_detected_umi.txt",
                    "_read_length.report.txt",
                    "_detected_dT.txt",
                    # FASTQ
                    ".chimeric_reads.fq.gz",
                    ".chimeric_segments.fq.gz",
                    ".non_chimeric_reads.fq.gz",
                    # BAM
                    "_chimeric_sorted.bam"),
    header = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
  )

  if(VALIDATE_CHIMERIC != TRUE){patterns <- patterns[!patterns$pattern %in% c("_mapping_chimerics.log.txt", "_chimeric_banList.txt", "_in_silico_chimerics.log.txt"),]}
  if(UMI_TYPE == "NONE"){patterns <- patterns[!patterns$pattern %in% c("vsearch_.*_UMI*.txt", "_get_UMI.report.txt", "_detected_umi.txt"),]}
  if(STRAND_SPECIFIC == FALSE | UMI_TYPE %in% c("MONOMER","TRIMER")){patterns <- patterns[!patterns$pattern %in% c("vsearch_*_dT.detection.txt", "_detected_dT.txt"),]}

  # Function to merge
  mergeFiles <- function(in_pattern = NULL, out_pattern = NULL, header = NULL){

    # FASTQ Files
    if(grepl("fq.gz$", in_pattern)){

      output_path <- paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_FLASHseqONT_output/", OUTPUT_SUFFIX, out_pattern)
      if(file.exists(output_path)){invisible(file.remove(output_path))}
      input_files <- list.files(OUTPUT_DIR, pattern = paste0(OUTPUT_SUFFIX, "_\\d+", out_pattern), full.names = TRUE)

      system(intern = FALSE,
             command = paste0("cat ", paste0(input_files, collapse = " "), " > ", output_path))

      # BAM FILES
    } else if(grepl("_sorted.bam$", in_pattern)){

      if(VALIDATE_CHIMERIC == TRUE){
        output_path <- paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_FLASHseqONT_output/", OUTPUT_SUFFIX, "_chimeric_sorted.bam")
        if(file.exists(output_path)){invisible(file.remove(output_path))}
        input_files <- list.files(OUTPUT_DIR, pattern = "_sorted.bam$", full.names = TRUE)

        system(intern = FALSE,
               command = paste0("samtools merge ", output_path, " ", paste0(input_files, collapse = " ")))
      }

      # Other text files
    } else if(grepl(".txt$", in_pattern)){

      output_path <- paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_FLASHseqONT_output/", OUTPUT_SUFFIX, out_pattern)

      if(file.exists(output_path)){invisible(file.remove(output_path))}

      if(in_pattern == "vsearch_*BC1*.txt"){
        input_files <- list.files(OUTPUT_DIR, pattern = paste0("vsearch_", OUTPUT_SUFFIX, "_BC1_\\d+", ".txt"), full.names = TRUE, recursive = TRUE)
      } else if(in_pattern == "vsearch_.*_UMI*.txt"){
        input_files <- list.files(OUTPUT_DIR, pattern = paste0("vsearch_", OUTPUT_SUFFIX, "_\\d+", "_UMI.txt"), full.names = TRUE, recursive = TRUE)
      } else if(in_pattern == "vsearch_*_dT.detection.txt"){
        input_files <- list.files(OUTPUT_DIR, pattern = paste0("vsearch_", OUTPUT_SUFFIX, "_\\d+", "_dT.detection.txt"), full.names = TRUE, recursive = TRUE)
      } else {
        input_files <- list.files(OUTPUT_DIR, pattern = paste0(OUTPUT_SUFFIX, "_\\d+", out_pattern), full.names = TRUE, recursive = TRUE)
      }

      if(length(input_files > 0)){

        if(header == TRUE){
          system(intern = FALSE,
                 command = paste0(
                   "head -n 1 ", input_files[1], " > ", output_path))

          system(intern = FALSE,
                 command = paste0(
                   "tail -n +2 -q ", paste0(input_files, collapse = " "), " >> ", output_path))
        } else {
          system(intern = FALSE,
                 command = paste0(
                   "cat ", paste0(input_files, collapse = " "), " > ", output_path))
        }

      } else {
        stop(paste0("STOP - Something is wrong, missing file for pattern: ", in_pattern))
      }

    }



  }

  # Run in parallel
  parallel::mclapply(1:nrow(patterns), mc.silent = silent, mc.cores = n_cores, function(i)
    mergeFiles(in_pattern = patterns$pattern[i], out_pattern = patterns$file_suffix[i], header = patterns$header[i])
  )


  # FASTQ - demultiplexed
  if(!silent){message("Demultiplexed FASTQ files")}

  demultiplexed_fastq <- list.files(paste0(OUTPUT_DIR, "/"), pattern = paste0(OUTPUT_SUFFIX, "_.*_demultiplexing_pass"), include.dirs = TRUE, full.names = TRUE, recursive = TRUE)
  fq.files <- as.vector(unlist(sapply(demultiplexed_fastq, function(x) list.files(x, full.names = TRUE, recursive = TRUE))))
  fq.files <- fq.files[grepl("_demultiplexed_.*\\.fq\\.gz", fq.files)]

  # Get the barcodes
  fq.bc <- str_extract(fq.files, pattern = "BC\\d+(?=\\.fq\\.gz$)|Chimeric|Undetermined")
  BC <- unique(unlist(fq.bc))

  # Write files
  output_path <- paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_FLASHseqONT_output/", OUTPUT_SUFFIX, "_demultiplexing_pass")
  if(file.exists(output_path)){invisible(unlink(output_path, recursive = TRUE))}
  invisible(dir.create(output_path))
  n_reads_written <- parallel::mclapply(mc.cores = THREADS, mc.silent = silent, BC, function(x)
    system(intern = FALSE,
           command = paste0("echo BC: ", x, " && cat ", paste0(fq.files[grepl(pattern = paste0(x, ".fq.gz$"), fq.files)], collapse = " "), " > ", paste0(output_path, "/", OUTPUT_SUFFIX, "_", x, "_demultiplexing_pass.fq.gz")))
  )

}


