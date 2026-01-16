#' @title Extract and Correct UMI sequences from reads (trimers)
#'
#' @description
#' Detect the UMI (trimer) position in the read (+/- 70%). Correct for indels / mutations - return the UMI associated to each read (+/- 85% efficacy)
#'
#' @param FASTQ character. "/path/to/the/input/fastq.gz". Gziped (expect '.gz') or not.
#' @param OUTPUT_DIR character. "/path/to/output/"
#' @param OUTPUT_SUFFIX character. sample id.
#' @param THREADS numeric. Number of threads. Will remove 1 for safety.
#' @param TRIMERS character vector.  c("TTC", "GGT", "AAA", "CCG").
#' @param WINDOW numeric. Detection windows for the UMI at the read start/ends.
#' @param dT_SEQ character. Should be set to "AAGCAGTGGTATCAACGCAGAGTNNNNNNNNNNNNNNNNNNNNNNNNCTTTTTTTTTT".
#' @param WHITELIST "ALL" or "CONTENT_AWARE". If ALL, use all the trimer combination for the whitelist. If CONTENT_AWARE, use all the unique perfect trimer combinations detected during the first pass.
#' @param silent boolean. If TRUE will print various messages.
#'
#' @export
get_and_correct_UMI_trimer <- function(FASTQ = NULL,
                                       OUTPUT_DIR = NULL,
                                       OUTPUT_SUFFIX = NULL,
                                       TRIMERS = c("TTC", "GGT", "AAA", "CCG"),
                                       dT_SEQ = "AAGCAGTGGTATCAACGCAGAGTNNNNNNNNNNNNNNNNNNNNNNNNCTTTTTTTTTT",
                                       THREADS = 1,
                                       WINDOW = 200,
                                       WHITELIST = "CONTENT_AWARE",
                                       silent = FALSE){

  # 1. Identify the dT position
  run_vsearch(FASTQ = FASTQ,
              OUTPUT_DIR = OUTPUT_DIR,
              FORWARD_SEQ = dT_SEQ,
              REVERSE_SEQ = dT_SEQ,
              SEQ_IDENTITY = 0.85,
              OUTPUT_SUFFIX = paste0(OUTPUT_SUFFIX, "_UMI"),
              MIN_SEQ_LENGTH = 30,
              WINDOW = WINDOW,
              THREADS = THREADS,
              silent = silent)

  # 2. Load and cleanup vsearch results
  vsearch_res <- suppressWarnings(read_tsv(paste0(OUTPUT_DIR, "/", "vsearch_", OUTPUT_SUFFIX, "_UMI.txt"), col_names = TRUE, show_col_types = FALSE)) %>%
    # Vsearch can find 2 hits (fwd / rev), select the most likely if there are two
    group_by(read_id) %>%
    dplyr::filter(bc_match_l == max(bc_match_l)) %>%
    dplyr::filter(identity == max(identity)) %>%
    dplyr::filter(n() == 1) %>%
    # Exclude missing bc
    dplyr::filter(barcode != "*") %>%
    ungroup() %>%
    select(barcode, read_id, query_seq_aln, bc_seq_aln) %>%
    # Same orientation for all sequences
    mutate(bc_seq_aln = ifelse(barcode == "BC1_rev", as.character(ShortRead::reverseComplement(Biostrings::DNAStringSet(bc_seq_aln))), bc_seq_aln),
           query_seq_aln = ifelse(barcode == "BC1_rev", as.character(ShortRead::reverseComplement(Biostrings::DNAStringSet(query_seq_aln))), query_seq_aln))

  # 3. Extract the UMI position
  bc_length <- nchar(str_extract(dT_SEQ, "N*N"))
  bc_positions <- str_locate_all(pattern = paste0(paste0(rep("N", bc_length), collapse = ""), "|", tolower(paste0(rep("N", bc_length), collapse = ""))),
                                 vsearch_res$bc_seq_aln)

  bc_sequences <- mapply(function(x,y) str_sub(string = y, start = x[,1], end = x[,2]), bc_positions, vsearch_res$query_seq_aln)
  bc_sequences[sapply(bc_sequences, length) == 0] <- NA_character_

  names(bc_sequences) <- vsearch_res$read_id

  # Tidy
  bc_sequences <- enframe(bc_sequences) %>%
    unnest(cols = c(value), keep_empty = TRUE) %>%
    mutate(orientation = vsearch_res$barcode,
           read_sequence = vsearch_res$query_seq_aln)
  colnames(bc_sequences) <- c("read_id", "value", "orientation", "read_sequence")

  # 4. Find the best reading frame

  # 4.1 FUNCTION: Split a DNA sequence by trimer. If 1 mismatch to the whitelist, return the closest else return the unmodified one.
  correct_trimer <- function(umi, trimers, THREADS){
    splitTrimer <- strsplit(umi, "(?<=.{3})", perl = TRUE)[[1]]
    dist_x <- stringdist::stringdistmatrix(splitTrimer, trimers, method = "h", nthread = THREADS)

    corrected <- sapply(1:nrow(dist_x), function(x) ifelse(any(dist_x[x,] <= 1), trimers[dist_x[x,] <= 1], splitTrimer[x]))

    return(corrected)
  }

  # 4.2 For each reading frame, return the correct (0-1 mimsatch) UMI and the number of trimer perfectly matching after correction
  trimer_frame <- list()
  for(i in 1:3){
    trimer_split <- lapply(bc_sequences$value, function(x) correct_trimer(substring(x,  1+i-1, nchar(x)), TRIMERS, THREADS))

    trimer_frame[[i]] <- data.frame(trimers = sapply(trimer_split, function(x) paste0(x,collapse = "")),
                                    unlist(sapply(trimer_split, function(x) sum(x %in% TRIMERS))))

    colnames(trimer_frame[[i]]) <- paste0(c("trimer_seq_", "nTrimer_"), i-1)

  }

  # 4.3. Tidy
  frames <- bind_cols(trimer_frame)
  frames$max_trimers <- apply(select(frames, starts_with("nTrimer_")), 1, function(x) max(x))
  frames$best_trimer <- apply(frames, 1, function(x) x[paste0("trimer_seq_", which.max(x[paste0("nTrimer_", 0:2)])-1)])

  bc_sequences <- bind_cols(bc_sequences, frames)

  # 5. Focus on the UMI that could not be corrected
  # Could contain some internal indels
  corrected_trimers <- dplyr::filter(bc_sequences, max_trimers == 8) %>%
    dplyr::select(read_id, best_trimer, orientation) %>%
    dplyr::rename(UMI = best_trimer, strand = orientation) %>%
    dplyr::mutate(confidence = "first_pass",
                  strand = ifelse(strand == "BC1_fwd", "+", "-"))

  UMI_toCorrect <- filter(bc_sequences, max_trimers != 8 & str_count(value, "-") < 2)

  # 5.1. Create a whitelist containing all possible UMI trimer combination
  if(WHITELIST == "ALL"){
    whitelist <- apply(expand_grid(TRIMERS, TRIMERS, TRIMERS, TRIMERS, TRIMERS, TRIMERS, TRIMERS, TRIMERS), 1, function(x) paste0(x, collapse = ""))
  } else if(WHITELIST == "CONTENT_AWARE"){
    whitelist <- unique(corrected_trimers$UMI)
  }

  names(whitelist) <- whitelist
  pUMI <- str_replace_all(UMI_toCorrect$best_trimer, "-", "N")
  names(pUMI) <- UMI_toCorrect$read_id
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(whitelist), paste0(OUTPUT_DIR, "/", "vsearch_tmp_", OUTPUT_SUFFIX, "_whitelist.fa"))
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(pUMI), paste0(OUTPUT_DIR, "/", "vsearch_tmp_", OUTPUT_SUFFIX, "_umi.fa"))

  # 5.2. RUN VSEARCH
  system(intern = FALSE,
         command = paste(
           "vsearch --usearch_global", paste0(OUTPUT_DIR, "/", "vsearch_tmp_", OUTPUT_SUFFIX, "_umi.fa"),
           "--threads", THREADS,
           "--minseqlength 24",
           "--maxaccepts 6",
           "--strand plus",
           "--wordlength 6",
           "--minwordmatches 9",
           "--userfields 'query+target+id+alnlen+mism+opens+qilo+qihi+qstrand+tilo+tihi+ql+tl+qrow+trow'",
           "--mincols 21",
           "--userout", paste0(OUTPUT_DIR, "/", "vsearch_tmp_", OUTPUT_SUFFIX, "_UMI.txt"),
           "--db", paste0(OUTPUT_DIR, "/", "vsearch_tmp_", OUTPUT_SUFFIX, "_whitelist.fa"),
           "--id", 0.85,
           ifelse(silent, "--quiet", "")
         )
  )

  system(intern = FALSE,
         command = paste(
           "sed -i",
           "'1i read_id	barcode	identity	bc_match_l	mism	gap_open	read_start	read_end	read_strand	bc_start	bc_end	query_l	bc_l	query_seq_aln	bc_seq_aln'",
           paste0(OUTPUT_DIR, "/", "vsearch_tmp_", OUTPUT_SUFFIX, "_UMI.txt")
         )
  )

  # 5.4. Mask ambiguous trimers from UMI with >85% match to the whitelist (= allow 2 mismatches per UMI)
  mask_diff <- function(bc) {

    diffs <- sapply(bc, function(x) strsplit(x, "(?<=.{3})", perl = TRUE)[[1]])

    index <- apply(diffs, 1, function(x) all(x == x[1]))

    return(paste0(ifelse(index, diffs[,1], "NNN"), collapse = ""))

  }

  vsearch_UMI_toCorrect <- read_tsv(paste0(OUTPUT_DIR, "/", "vsearch_tmp_", OUTPUT_SUFFIX, "_UMI.txt"), col_names = TRUE, show_col_types = FALSE)

  if(nrow(vsearch_UMI_toCorrect) > 0){

    vsearch_UMI_toCorrect <- vsearch_UMI_toCorrect %>%
      group_by(read_id) %>%
      # Select end-to-end matches
      filter(bc_match_l == max(bc_match_l)) %>%
      # If multiple possibilities, select highest identity scores
      filter(identity == max(identity)) %>%
      # Compute number of match per UMI
      mutate(n = n()) %>%
      # If more than one good hit is reported for a UMI - compare them and mask their differences with Ns
      # For instance: "CCGGGTTTCGGTTTCAAATTCTTC" "CCGGGTTTCGGTAAAAAATTCTTC"
      # Becomes: "CCGGGTTTCGGTNNNAAATTCTTC"
      # Also works to correct match where the first/last trimer is not reported by VSEARCH - to complete the UMI to 24bp
      reframe(UMI = ifelse(n != 1, mask_diff(barcode), barcode)) %>%
      ungroup() %>%
      # Remove duplicates
      distinct() %>%
      # Filter out elements with too many NNN
      filter(str_count(pattern = "NNN", UMI) <3) %>%
      mutate(confidence = "second_pass") %>%
      # Add strand orientation
      left_join(select(bc_sequences, read_id, orientation), by = "read_id") %>%
      dplyr::rename(strand = orientation) %>%
      mutate(strand = ifelse(strand == "BC1_fwd", "+", "-"))
  } else {
    message("No match found")

  }

  # 6. Some stats
  if(!silent){
    message(paste0("n_reads UMI total (after first pass): ", nrow(bc_sequences) - nrow(UMI_toCorrect), " (", floor(100*(nrow(bc_sequences) - nrow(UMI_toCorrect)) / nrow(bc_sequences)), "%)"))
    message(paste0("n_reads UMI total (after second pass): ", nrow(bc_sequences) - nrow(UMI_toCorrect) + nrow(vsearch_UMI_toCorrect), " (", floor(100*(nrow(vsearch_UMI_toCorrect) + nrow(bc_sequences) - nrow(UMI_toCorrect)) / nrow(bc_sequences)), "%)"))
  }
  # 7. Return results
  umi_list <- bind_rows(corrected_trimers,
                        vsearch_UMI_toCorrect) %>%
  # 7.1 Collapse trimers for umi_tools
  mutate(umi_collapsed = sapply(UMI, function(x) paste0(strsplit(x, "")[[1]][seq(1, 24, by = 3)], collapse = "")))

  # 8. Cleanup
  unlink(paste0(OUTPUT_DIR, "/", "vsearch_tmp_", OUTPUT_SUFFIX, "*.fa"))
  unlink(paste0(OUTPUT_DIR, "/", "vsearch_tmp_", OUTPUT_SUFFIX, "_UMI.txt"))

  # 9. Write Report
  write.table(file = paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_get_UMI.report.txt"), sep="\t", row.names = FALSE, quote = FALSE,
              data.frame(
                ID = OUTPUT_SUFFIX,
                n_vsearch_hit_umidTseq = length(unique(vsearch_res$read_id)),
                n_umi_first_pass = length(unique(corrected_trimers$read_id)),
                n_umi_first_pass_fwdstrand = sum(corrected_trimers$strand == "+"),
                n_umi_second_pass = length(unique(vsearch_UMI_toCorrect$read_id)),
                n_umi_second_pass_fwdstrand = sum(vsearch_UMI_toCorrect$strand == "+"),
                n_umi_masked_bases = sum(str_detect(pattern = "NNN", umi_list$UMI))) %>%
                mutate(
                  p_umi_first_pass = round(100*n_umi_first_pass/n_vsearch_hit_umidTseq, 2),
                  p_umi_second_pass = round(100*n_umi_second_pass/n_vsearch_hit_umidTseq, 2),
                  p_umi_total = round(100*(p_umi_first_pass+p_umi_second_pass)/n_vsearch_hit_umidTseq, 2),
                  p_umi_masked_bases = round(100*n_umi_masked_bases/n_vsearch_hit_umidTseq, 2),
                  correct_trimer_0 = sum(bc_sequences$max_trimers == 0),
                  correct_trimer_1 = sum(bc_sequences$max_trimers == 1),
                  correct_trimer_2 = sum(bc_sequences$max_trimers == 2),
                  correct_trimer_3 = sum(bc_sequences$max_trimers == 3),
                  correct_trimer_4 = sum(bc_sequences$max_trimers == 4),
                  correct_trimer_5 = sum(bc_sequences$max_trimers == 5),
                  correct_trimer_6 = sum(bc_sequences$max_trimers == 6),
                  correct_trimer_7 = sum(bc_sequences$max_trimers == 7),
                  correct_trimer_8 = sum(bc_sequences$max_trimers == 8)
                )
  )

  # 10. Write UMI
  write.table(umi_list, paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_detected_umi.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

  # return(umi_list)
}

#' @title Extract and Correct UMI sequences from reads
#'
#' @description
#' Detect the UMI (trimer) position in the read.
#'
#' NB: Do not modify the minimal seq identity - it was chosen to maximise the recovery of UMI while excluding oligo-dT without UMI.
#' Will report <0.2% of matching reads in samples without oligo-dT-UMI.
#'
#' @param FASTQ character. "/path/to/the/input/fastq.gz". Gziped (expect '.gz') or not.
#' @param OUTPUT_DIR character. "/path/to/output/"
#' @param OUTPUT_SUFFIX character. sample id.
#' @param THREADS numeric. Number of threads. Will remove 1 for safety.
#' @param dT_SEQ character. Should be set to "AAGCAGTGGTATCAACGCAGAGTCTTTNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTT".
#' @param MAX_GAP_LENGTH numeric. Max gap length in detected dT_SEQ.
#' @param MAX_GAPS_N numeric. Max number of gaps allowed in detected dT_SEQ.
#' @param MAX_MISMATCH_N numeric. Max number of mismatches allowed in detected dT_SEQ.
#' @param silent boolean. If TRUE will print various messages.
#'
#' @export
get_UMI <- function(FASTQ = NULL,
                    OUTPUT_DIR = NULL,
                    OUTPUT_SUFFIX = NULL,
                    dT_SEQ = "AAGCAGTGGTATCAACGCAGAGTACNNNNNNNNATACTGACGCTTT",
                    THREADS = 1,
                    MAX_GAP_LENGTH = 2,
                    MAX_GAPS_N = 2,
                    WINDOW = 200,
                    MAX_MISMATCH_N = 6,
                    silent = FALSE){

  # 1. Identify the dT position
  run_vsearch(FASTQ = FASTQ,
              OUTPUT_DIR = OUTPUT_DIR,
              FORWARD_SEQ = dT_SEQ,
              REVERSE_SEQ = dT_SEQ,
              SEQ_IDENTITY = 0.86,
              OUTPUT_SUFFIX = paste0(OUTPUT_SUFFIX, "_UMI"),
              MIN_SEQ_LENGTH = 30,
              WINDOW = WINDOW,
              THREADS = THREADS,
              silent = silent)

  # 2. Load and cleanup vsearch results
  vsearch_res <- suppressWarnings(vroom::vroom(paste0(OUTPUT_DIR, "/vsearch_", OUTPUT_SUFFIX, "_UMI.txt"), delim = "\t", col_names = TRUE, show_col_types = FALSE))

  if(nrow(vsearch_res) > 10){

    gap_length_max <- function(x = NULL){
      max((rle(strsplit(x, "")[[1]])$lengths)[(rle(strsplit(x, "")[[1]])$value == "-")]) %>% return
    }

    vsearch_res <- vsearch_res %>%
      group_by(read_id) %>%
      dplyr::filter(gap_open <= MAX_GAPS_N,
                    mism <= MAX_MISMATCH_N,
                    # Exclude match were identity is good but length is insufficient (= could be lone ISPCR)
                    bc_match_l > floor(nchar(dT_SEQ) * 0.85),
                    # Exclude missing bc if any (should not)
                    barcode != "*") %>%
      # Vsearch can find 2 hits (fwd / rev), in that case, remove both
      dplyr::filter(n() == 1) %>%
      ungroup() %>%
      select(barcode, read_id, query_seq_aln, bc_seq_aln) %>%
      # Same orientation for all sequences
      mutate(bc_seq_aln = ifelse(barcode == "BC1_rev", as.character(ShortRead::reverseComplement(Biostrings::DNAStringSet(bc_seq_aln))), bc_seq_aln),
             query_seq_aln = ifelse(barcode == "BC1_rev", as.character(ShortRead::reverseComplement(Biostrings::DNAStringSet(query_seq_aln))), query_seq_aln))

    # Exclude long gaps
    # Filter out reads with large / numerous gaps
    vsearch_res$max_gap_length <- suppressWarnings(sapply(vsearch_res$query_seq_aln, function(x) gap_length_max(x)))
    vsearch_res$max_gap_length <- ifelse(vsearch_res$max_gap_length == -Inf, 0, vsearch_res$max_gap_length)

    vsearch_res <- filter(vsearch_res, max_gap_length <= MAX_GAP_LENGTH)

    # 3. Extract the UMI position
    bc_length <- nchar(str_extract(dT_SEQ, "N*N"))
    bc_positions <- str_locate_all(pattern = paste0(paste0(rep("N", bc_length), collapse = ""), "|", tolower(paste0(rep("N", bc_length), collapse = ""))),
                                   vsearch_res$bc_seq_aln)

    bc_sequences <- mapply(function(x,y) str_sub(string = y, start = x[,1], end = x[,2]), bc_positions, vsearch_res$query_seq_aln)
    bc_sequences[sapply(bc_sequences, length) == 0] <- NA_character_

    names(bc_sequences) <- vsearch_res$read_id

    # Tidy
    bc_sequences <- enframe(bc_sequences) %>%
      unnest(cols = c(value), keep_empty = TRUE) %>%
      mutate(orientation = vsearch_res$barcode,
             read_sequence = vsearch_res$query_seq_aln)
    colnames(bc_sequences) <- c("read_id", "value", "orientation", "read_sequence")

    # Mask DEL
    bc_sequences$value <- str_replace_all(bc_sequences$value, "-", "N")

    # Return
    umi_list <- bc_sequences %>%
      mutate(strand = ifelse(orientation == "BC1_fwd", "+", "-")) %>%
      select(read_id, value, strand) %>%
      dplyr::rename("UMI" = value)

    # Write Report
    #message(paste0("Write: ", paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_get_UMI.report.txt")))
    write.table(file = paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_get_UMI.report.txt"), sep = "\t", row.names = FALSE, quote = FALSE,
                data.frame(
                  ID = OUTPUT_SUFFIX,
                  n_vsearch_hit_umidTseq = length(unique(vsearch_res$read_id)),
                  n_umi = length(unique(umi_list$read_id)),
                  n_umi_fwdstrand = sum(umi_list$strand == "+"),
                  n_umi_masked_bases = sum(str_detect(pattern = "N", umi_list$UMI))) %>%
                  mutate(
                    p_umi = round(100*n_umi/n_vsearch_hit_umidTseq, 2),
                    p_umi_masked_bases = round(100*n_umi_masked_bases/n_vsearch_hit_umidTseq, 2)
                  )

    )

    # 10. Write UMI
    #message(paste0("Write: ", paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_detected_umi.txt")))
    #message(paste0("Found UMI:", umi_list), " of ", nrow(umi_list))
    write.table(umi_list, paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_detected_umi.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

    # return(umi_list)

  } else {
    message("Not enough UMI reads detected")
    # return(tibble())
  }

  # Clean-up
  unlink(paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "*.fa"))
}



