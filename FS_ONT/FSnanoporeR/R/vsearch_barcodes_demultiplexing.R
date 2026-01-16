#' filter Trim and demultiplex reads
#'
#' THIS FUNCTION WILL LIKELY BE SPLIT IN MULTIPLE PARTS - BECAME TOO BIG.
#'
#' @description
#' Take the output from 'extractBarcodes()' and remove the chimeric indexes combinations and optionally select specific barcodes from the whitelist.
#' Create individual files containing the read IDs associated to a specific barcode
#'
#' @param BARCODES_PATH data.frame. Output from 'extractBarcodes()'.
#'
#' @param FASTQ character. "/path/to/the/input/fastq.gz". Gziped (.gz) or not.
#' @param THREADS numeric. Number of threads to run bbmap.
#' @param OUTPUT_DIR character. /path/to/output_dir/.
#' @param OUTPUT_SUFFIX character. Sample id.
#' @param silent boolean. Display messages or not.
#' @param N_READS numeric. For the demultiplexing, iterate through chunks of N_READS (see ShortRead::FastqStreamer).
#'
#' @param UMI_TYPE, character. Accepts "NONE", "MONOMER", "TRIMER" or default to no UMI.
#' @param TRIM boolean. (Beta) If TRUE, will search for the PRIMER_SEQ and trim anything from it to the read ends. Currently the trimming will leave CCG (TSO) or T(30) (oligo-dT) at the end of the reads.
#' @param VSEARCH_BARCODE_RESULTS character. If trimming, get the barcode positions from vsearch.
#' @param BANLIST character. "/path/to/banlist_dataframe.txt". Created with 'chimeric_mapping_validation()'
#' @param WINDOW numeric. Run vsearch only on the read start/end windows. Must be the same as WINDOW in vsearch_run.
#'
#' @import tidyverse
#' @importFrom ShortRead yield
#' @importFrom ShortRead FastqStreamer
#' @importFrom ShortRead writeFastq
#' @importFrom ShortRead readFastq
#' @import parallel
#'
#' @returns Write demultiplexed FASTQ files to OUTPUT_DIR
#'
#' @export
filter_trim_demultiplex_OLD <- function(BARCODES_PATH = NULL,
                                    FASTQ = NULL,
                                    OUTPUT_DIR = NULL,
                                    OUTPUT_SUFFIX = NULL,
                                    THREADS = 1,
                                    N_READS = 10000,
                                    UMI_TYPE = "NONE",
                                    STRAND_SPECIFIC = FALSE,
                                    WINDOW = NULL,
                                    BANLIST = NULL,
                                    TRIM = TRUE,
                                    VSEARCH_BARCODE_RESULTS = NULL,
                                    silent = FALSE){

  if(!dir.exists(OUTPUT_DIR)){dir.create(OUTPUT_DIR)}

  # 1. Remove chimeric index combinations
  if(!silent){message("1. Remove chimeric barcode combinations")}

  bc_filtered <- read_tsv(BARCODES_PATH, show_col_types = FALSE) %>%
    select(read_id, index_unique)

  # 3. Optional - Filter out banned reads (chimeric)
  # For now ban any chimeric that have no match with the read mapping or in-silico chimerics
  if(file.exists(BANLIST)){
    banlist <- read_tsv(BANLIST, show_col_types = FALSE)
    bc_filtered <- filter(bc_filtered, !read_id %in% banlist$read_id)
  }

  # 5. Create folders
  files_to_create <- paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_demultiplexed_", unique(bc_filtered$index_unique), ".fq.gz")
  if(any(file.exists(files_to_create))){
    if(!silent){message("Erase previously created FASTQ to avoid appending the new results to them")}
    invisible(file.remove(files_to_create))
  }

  # 6. Load the FASTQ file
  if(!silent){message(paste0("3. Load FASTQ file - SLOW "))}
  fq.opened <- ShortRead::readFastq(FASTQ)

  total_reads <- length(fq.opened)

  # 7. UMI detection and read trimming

  # 7.1. Detect UMI
  if(UMI_TYPE %in% c("TRIMER", "MONOMER")){
    if(!silent){message(paste0("4. Get UMIs and add barcode information"))}
    # 7.1.1. Regular UMI
    if(UMI_TYPE == "MONOMER"){
      message("Monomeric UMI Selected")
      umi_list <- get_UMI(FASTQ = FASTQ,
                          OUTPUT_DIR = OUTPUT_DIR,
                          OUTPUT_SUFFIX = OUTPUT_SUFFIX,
                          WIN = WINDOW,
                          dT_SEQ = "AAGCAGTGGTATCAACGCAGAGTACNNNNNNNNATACTGACGCTTT",
                          THREADS = THREADS)

      # 7.1.2. Trimer UMI (see Sun et al 2023)
    } else if(UMI_TYPE == "TRIMER"){
      umi_list <- get_and_correct_UMI_trimer(FASTQ = FASTQ,
                                             OUTPUT_DIR = OUTPUT_DIR,
                                             OUTPUT_SUFFIX = OUTPUT_SUFFIX,
                                             TRIMERS = c("TTC", "GGT", "AAA", "CCG"),
                                             dT_SEQ = "AAGCAGTGGTATCAACGCAGAGTNNNNNNNNNNNNNNNNNNNNNNNNCTTTTTTTTTTTTT",
                                             WIN = WINDOW,
                                             THREADS = 1,
                                             silent = FALSE)
    }


    # 7.2. Append UMI to read ID (UMI-tools)
    # Ideally would have it in unmapped bam but not accepted by minimap2
    fq.opened <- tibble(read_id = str_split(as.character(ShortRead::id(fq.opened)), " ", simplify = TRUE)[,1],
                        read_seq = as.character(ShortRead::sread(fq.opened)),
                        read_qual = as.character(fq.opened@quality@quality),
                        read_length = ShortRead::width(fq.opened)) %>%
      mutate(read_id = str_split(read_id, " ", simplify = TRUE)[,1]) %>%
      # Add barcode information
      left_join(bc_filtered, by = "read_id") %>%
      # Add UMI information
      left_join(umi_list, by ="read_id") %>%
      dplyr::rename("orientation" = "strand") %>%
      mutate(orientation = ifelse(is.na(orientation), "*", orientation))

    if(UMI_TYPE == "MONOMER"){
      umi.seq <- fq.opened$UMI
    }else{
      umi.seq <- fq.opened$umi_collapsed
    }

    # 7.3. Update read IDs with UMI info
    fq.opened <- fq.opened %>%
      # remove cases where there are no UMI associated
      mutate(read_id.updated = paste0(read_id, "_", umi.seq),
             read_id.updated = str_replace(read_id.updated, "_NA", "")) %>%
      # Select relevant columns
      dplyr::select(read_id, read_seq, read_qual, read_length, index_unique, read_id.updated, orientation)

    umi_list <- NULL

  } else {
    # If no UMI are present
    if(!silent){message(paste0("4. Tidy and add barcode information"))}
    fq.opened <- tibble(read_id = str_split(as.character(ShortRead::id(fq.opened)), " ", simplify = TRUE)[,1],
                        read_seq = as.character(ShortRead::sread(fq.opened)),
                        read_qual = as.character(fq.opened@quality@quality),
                        read_length = ShortRead::width(fq.opened)) %>%
      # Add barcode information
      left_join(bc_filtered, by = "read_id") %>%
      # Select relevant columns
      dplyr::select(read_id, read_seq, read_qual, read_length, index_unique)

    if(STRAND_SPECIFIC == TRUE){
      # Add the read orientation based on the oligo-dT presence
      orientation <- determine_strand_orientation(FASTQ = FASTQ,
                                                  OUTPUT_DIR = OUTPUT_DIR,
                                                  OUTPUT_SUFFIX = OUTPUT_SUFFIX,
                                                  dT_SEQ = "AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTT",
                                                  THREADS = 1,
                                                  WINDOW = WINDOW,
                                                  silent = silent)

      fq.opened <- left_join(fq.opened, select(orientation, -barcode), by = "read_id") %>%
        mutate(orientation = ifelse(is.na(orientation), "*", orientation))
    }

  }

  # 8. Remove reads without barcodes
  if(!silent){message(paste0("5. Filter out reads without barcodes"))}
  fq.opened <- filter(fq.opened, !is.na(index_unique))

  # 9. Read Trimming
  if(TRIM){
    if(!silent){message(paste0("6. Optional - Trim sequences"))}

    # 9.1. Load the barcode position information
    vsearch_bc_res <- suppressWarnings(read_tsv(VSEARCH_BARCODE_RESULTS, show_col_types = FALSE)) %>%
      filter(read_id %in% fq.opened$read_id)

    if(UMI_TYPE %in% c("TRIMER", "MONOMER")){
      # 9.2.1. Load - if available - the UMI position information
      vsearch_umi_res <- suppressWarnings(read_tsv(paste0(OUTPUT_DIR, "/", "vsearch_", OUTPUT_SUFFIX, "_UMI.txt"), show_col_types = FALSE)) %>%
        filter(read_id %in% fq.opened$read_id)

      # 9.2.2. Merge both
      vsearch_bc_res <- bind_rows(vsearch_bc_res, vsearch_umi_res) %>%
        # Add Read Lengths
        left_join(dplyr::select(fq.opened, read_id, read_length),
                  by = "read_id")

      vsearch_bc_res_readStarts <- filter(vsearch_bc_res, read_section == "START") %>%
        group_by(read_id) %>%
        summarise(trim_pos = unique(max(read_end))) %>%
        ungroup()

      vsearch_bc_res_readEnds <- filter(vsearch_bc_res, read_section == "END") %>%
        group_by(read_id) %>%
        summarise(trim_pos = unique(as.numeric(read_length)-(query_l-min(read_start)+1))) %>%
        ungroup()

      vsearch_bc_res <- left_join(vsearch_bc_res_readEnds, vsearch_bc_res_readStarts, by = "read_id") %>%
        dplyr::rename("trim_end_read" = trim_pos.x,  "trim_start_read" = trim_pos.y)

      vsearch_bc_res_readEnds <- vsearch_bc_res_readStarts <- vsearch_umi_res <- NULL

    } else {
      # Get barcode start/end
      vsearch_bc_res <- vsearch_bc_res %>%
        # Add Read Lengths
        left_join(dplyr::select(fq.opened, read_id, read_length),
                  by = "read_id") %>%
        group_by(read_id, read_section) %>%
        # Adjust trimming positions
        summarise(
          trim_pos = ifelse(read_section == "START", read_end, read_length-(query_l-read_start+1))
        ) %>%
        ungroup() %>%
        # Open df
        pivot_wider(names_from = read_section, values_from = trim_pos) %>%
        dplyr::rename("trim_end_read" = END,  "trim_start_read" = START)
    }

    fq.opened <- left_join(fq.opened,
                           vsearch_bc_res, by = c("read_id")) %>%
      # Correct for missing barcode information
      mutate(
        trim_start_read = ifelse(is.na(trim_start_read), 0, trim_start_read),
        trim_end_read = ifelse(is.na(trim_end_read), read_length, trim_end_read)
      ) %>%
      # Trim Reads
      mutate(
        read_seq = str_sub(read_seq, start = trim_start_read, end = trim_end_read),
        read_qual = str_sub(read_qual, start = trim_start_read, end = trim_end_read)
      )

    vsearch_bc_res <- NULL

  }

  # 10. Write separated FASTQ
  if(!silent){message(paste0("7. Demultiplexing"))}
  unique_bc <- unique(fq.opened$index_unique)

  toFASTQ <- function(FQ = NULL, BC = NULL, STRAND_SPECIFIC = FALSE){

    fq.opened.tmp <- filter(fq.opened, index_unique == BC)
    # Reverse complement
    if(STRAND_SPECIFIC == TRUE){
      fq.opened.tmp <- mutate(fq.opened.tmp,
                              read_seq = ifelse(orientation == "+", Biostrings::reverseComplement(Biostrings::DNAStringSet(read_seq)), read_seq))
    }

    fq <- ShortRead::ShortReadQ(
      sread = Biostrings::DNAStringSet(fq.opened.tmp$read_seq),
      id = Biostrings::BStringSet(fq.opened.tmp$read_id),
      quality = Biostrings::BStringSet(fq.opened.tmp$read_qual)
    )

    return(fq)
  }

  # Use the read ids with UMI attached if necessary
  if(UMI_TYPE %in% c("MONOMER", "TRIMER")){
    fq.opened$read_id <- fq.opened$read_id.updated
  }

  unlink(paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_demultiplexed_*.fq.gz"))
  n_reads_written <- parallel::mclapply(mc.cores = THREADS, mc.silent = TRUE, unique_bc, function(x)
    ShortRead::writeFastq(object = toFASTQ(FQ = fq.opened, BC = x, STRAND_SPECIFIC = STRAND_SPECIFIC),
                          file = paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_demultiplexed_", x, ".fq.gz"),
                          mode = "w",
                          full=FALSE)
  )


  # 12. Reports
  if(!silent){message(paste0("8. Create summary report"))}
  write.table(file = paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_demultiplex.report.txt"), sep="\t", row.names = FALSE, quote = FALSE,
              fq.opened %>%
                group_by(index_unique) %>%
                summarise(n_reads_before_demultiplexing = n(),
                          median_width = median(as.numeric(read_length)),
                          sd_width = sd(as.numeric(read_length))
                ) %>%
                mutate(ID = OUTPUT_SUFFIX) %>%
                dplyr::rename("Barcode" = "index_unique")
  )

  write.table(file = paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_read_length.report.txt"), sep="\t", row.names = FALSE, quote = FALSE,
              data.frame(
                read_id = fq.opened$read_id,
                read_width = fq.opened$read_length,
                read_bp_score_total = ShortRead::alphabetScore(ShortRead::FastqQuality(Biostrings::BStringSet(fq.opened$read_qual))),
                cell_barcode = fq.opened$index_unique
              )
  )

  # 13. End
  if(!silent){message("Finished")}

}





#' determine_strand_orientation
#'
#' @description
#' Detect the position of the dT in the read. Useful to reorient the reads afterwards.
#'
#' @param FASTQ character. "/path/to/the/input/fastq.gz". Gziped (.gz) or not.
#' @param THREADS numeric. Number of threads to run bbmap.
#' @param OUTPUT_DIR character. /path/to/output_dir/.
#' @param OUTPUT_SUFFIX character. Sample id.
#' @param silent boolean. Display messages or not.
#' @param dT_SEQ character. Sequence of the dT (no UMI / Ns allowed)
#' @param WINDOW numeric. Up/down window in the reads to look for the dT
#'
#' @import tidyverse
#'
#' @returns *_detected_dT.txt object
#'
#' @export
determine_strand_orientation <- function(FASTQ = NULL,
                                         OUTPUT_DIR = NULL,
                                         OUTPUT_SUFFIX = NULL,
                                         dT_SEQ = "AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTT",
                                         THREADS = 1,
                                         WINDOW = 200,
                                         silent = FALSE){

  # 1. Identify the dT position
  run_vsearch(FASTQ = FASTQ,
                OUTPUT_DIR = OUTPUT_DIR,
                FORWARD_SEQ = dT_SEQ,
                REVERSE_SEQ = dT_SEQ,
                SEQ_IDENTITY = 0.5,
                OUTPUT_SUFFIX = paste0(OUTPUT_SUFFIX, "_dT.detection"),
                MIN_SEQ_LENGTH = 30,
                WINDOW = WINDOW,
                THREADS = THREADS,
                silent = silent)

  # 2. Filter results
  # V1 - works well but dT and TSO are a bit too close.
  # Idea: V2 dT should integrate a few extra nucleotides to help detection.
  dT.positions <- read_tsv(paste0(OUTPUT_DIR, "/vsearch_", OUTPUT_SUFFIX, "_dT.detection.txt"), show_col_types = FALSE) %>%
    rowwise() %>%
    mutate(polyT = str_count(tolower(query_seq_aln), "tttt|aaaa") >= 3 | str_count(tolower(query_seq_aln), "ttttt|aaaaa") >= 2) %>%
    # 0.26% of the reads were likely from TSO binding close to a poly-T or poly-A seq.
    filter(!grepl("AGTACGGG|CCCGTACT", toupper(query_seq_aln)),
           polyT == TRUE) %>%
    # Vsearch can find 2 hits (fwd / rev), select the most likely if there are two
    group_by(read_id) %>%
    dplyr::filter(bc_match_l == max(bc_match_l)) %>%
    dplyr::filter(identity == max(identity)) %>%
    dplyr::filter(n() == 1) %>%
    ungroup() %>%
    select(read_id, barcode) %>%
    mutate(orientation = ifelse(barcode == "BC1_fwd", "+", "-"))

  write.table(dT.positions, paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_detected_dT.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

  # return(dT.positions)

}
