#' Extract barcodes from VSEARCH results
#'
#' @description
#' Extract the barcode sequences from the VSEARCH results ('run_vsearch()').
#' Compare them to a index whitelist. Sequencing / PCR errors are corrected by taking the closest matching barcode sequence.
#' Discard reads with barcodes too distant from the whitelist (>MAX_LEVENSTEIN). Given the current barcode design, it is not recommended to push >=3.
#'
#' @param VSEARCH_RES_PATH character. "/path/to/the/input/vsearch_results.txt".
#' @param WHITELIST_PATH character. "/path/to/the/whitelist_index.fasta".
#' @param THREADS numeric. Thread number to calculate the distances.
#' @param SELECTED_BC_IDs character. Vector of the barcode IDs to consider. Ex. c("BC1", "BC2", ...).
#' @param MIN_READS numeric. Output only indexes supported with at least MIN_READS
#' @param OUTPUT_DIR character. "/path/to/output/".
#' @param PRIMER_SEQ character.
#' @param OUTPUT_SUFFIX character.
#' @param silent boolean. Display messages / stats.
#'
#' @import tidyverse
#' @import ShortRead
#' @import stringdist
#'
#' @returns
#' dataframe associating the most likely index to each read.
#'
#' @export
extractBarcodes <- function(VSEARCH_RES_PATH = NULL,
                            silent = FALSE,
                            WHITELIST_PATH = NULL,
                            SELECTED_BC_IDs = NULL,
                            THREADS = 4,
                            MIN_READS = 100,
                            PRIMER_SEQ = NULL,
                            OUTPUT_DIR = NULL,
                            OUTPUT_SUFFIX = NULL){

  # 1. Load
  whitelist <- as.character(Biostrings::readDNAStringSet(WHITELIST_PATH))

  if(!is.null(SELECTED_BC_IDs)){
    if(!silent){message("SELECTED_BC ON - search only for specific barcodes")}
    whitelist <- whitelist[names(whitelist) %in% SELECTED_BC_IDs]
  }

  # NB: One warning about parsing issue for the file. Can be disregarded.
  vsearch_res <- suppressWarnings(read_tsv(VSEARCH_RES_PATH, col_names = TRUE, show_col_types = FALSE))

  # 2. Determine the maximal acceptable levenstein distance
  bc_distances <- stringdist::stringdistmatrix(whitelist, whitelist, nthread = 1, method = "lv")
  bc_distances <- tibble(
    bc_id = names(whitelist),
    bc_seq = whitelist,
    max_allowed_dist = apply(bc_distances, 2, function(x) min(x[x != 0]))-1
  )
  MAX_LEVENSTEIN <- min(bc_distances$max_allowed_dist)

  if(stringdist::stringdistmatrix(whitelist, nthread = 1, method = "lv") %>% min <= MAX_LEVENSTEIN){
    message("MAX_LEVENSTEIN is inferior or equal to the smallest levenstein distance of the indexes. To avoid collisions please decrease this value.")
  }

  # 3. Filter out reads without detected primer
  total_reads <- length(unique(vsearch_res$read_id))
  barcoded_reads <- length(unique(vsearch_res$read_id[vsearch_res$barcode != "*"]))
  if(!silent){message(paste0("Unique reads: ", total_reads))}

  vsearch_res <- dplyr::filter(vsearch_res, barcode != "*") %>%
    select(barcode, read_id, query_seq_aln, bc_seq_aln) %>%
    mutate(bc_seq_aln = ifelse(barcode == "BC1_rev", as.character(ShortRead::reverseComplement(Biostrings::DNAStringSet(bc_seq_aln))), bc_seq_aln),
           query_seq_aln = ifelse(barcode == "BC1_rev", as.character(ShortRead::reverseComplement(Biostrings::DNAStringSet(query_seq_aln))), query_seq_aln))

  # 4. Extract BC position / sequence
  # Take it into account for the distance
  bc_length <- nchar(str_extract(PRIMER_SEQ, "N*N"))
  bc_positions <- str_locate_all(pattern = paste0(paste0(rep("N", bc_length), collapse = ""), "|", tolower(paste0(rep("N", bc_length), collapse = ""))),
                                 vsearch_res$bc_seq_aln)

  bc_sequences <- mapply(function(x,y) str_sub(string = y, start = x[,1], end = x[,2]), bc_positions, vsearch_res$query_seq_aln)
  bc_sequences[sapply(bc_sequences, length) == 0] <- NA_character_

  names(bc_sequences) <- vsearch_res$read_id

  # 5. Tidy
  bc_sequences <- enframe(bc_sequences) %>%
    unnest(cols = c(value), keep_empty = TRUE) %>%
    mutate(orientation = vsearch_res$barcode,
           read_sequence = vsearch_res$query_seq_aln) %>%
    # Filter out when vsearch does not find the NNNN...NNNN barcode sequence
    filter(!is.na(value))
  colnames(bc_sequences) <- c("read_id", "value", "orientation", "read_sequence")

  # 6. Compare extracted sequences to the whitelist

  # 6.1. Separate perfect match from partial
  bc_sequences_perfect <- filter(bc_sequences, value %in% whitelist) %>%
    mutate(min_lv = 0,
           is_real_bc = TRUE,
           TYPE = "PERFECT") %>%
    left_join(as.data.frame(whitelist) %>% rownames_to_column("best_index"), by = c("value" = "whitelist"))

  bc_sequences_tocheck <- filter(bc_sequences,
                                 !value %in% whitelist,
                                 str_count(pattern = "-", value) <= 3)

  # 6.2. Get the levenstein distance of each detected BC to the whitelist
  lv_dist <- sapply(bc_sequences_tocheck$value, function(x)
    list(stringdist::stringdist(whitelist, x, method = "lv", nthread = THREADS)))

  # 6.3. Get the bc index with the smallest distance to an element of the whitelist
  bc_sequences_tocheck <- bc_sequences_tocheck %>%
    mutate(
      # Levenstein distance from one bc to the whitelist
      lv_dist = lv_dist,
      # Minimal LV distance
      min_lv = purrr::map(lv_dist, min),
      # How many close bc have the same min dist ?
      n_same_dist = purrr::map2_int(lv_dist, min_lv, ~sum(.x == .y)),
      # Based on the smallest distance what are the potential BC ?
      possible_indexes = purrr::map2(lv_dist, min_lv, ~names(whitelist)[.x == .y]),
      max_allowed_lv = MAX_LEVENSTEIN)

  # 6.4. Recombine the information of fwd/rev barcodes from the same read.
  real_bc <- bc_sequences_tocheck %>%
    unnest(cols = min_lv, keep_empty = TRUE) %>%
    # Filter out BC with too high lv dist and those matching >=2 barcodes from the whitelist
    filter(min_lv < max_allowed_lv & n_same_dist < 3) %>%
    select(-read_sequence, -lv_dist, -max_allowed_lv) %>%
    mutate(TYPE = "MISMATCHES") %>%
    select(-value, -min_lv) %>%
    dplyr::rename("best_index" = "possible_indexes") %>%
    bind_rows(bc_sequences_perfect %>%
                select(read_id, orientation, best_index, TYPE) %>%
                mutate(best_index = as.list(best_index),
                       n_same_dist = 1)
    ) %>%
    # Put in parallel the correction TYPE and indexes
    pivot_wider(
      id_cols = read_id,
      names_from = orientation,
      values_from = c("TYPE", "best_index", "n_same_dist"), names_sep = "-"
    ) %>%
    # Determine how many which and how many barcodes are in common between fwd/rev
    mutate(
      intersect_bc = purrr::map2(`best_index-BC1_fwd`, `best_index-BC1_rev`, ~intersect(.x, .y)),
      n_common = purrr::map(intersect_bc, ~length(.x))) %>%
    unnest(n_common, keep_empty = TRUE) %>%
    # Decision tree
    #
    mutate(
      best_index =
        case_when(
          # Perfect match & same bc both ends
          `TYPE-BC1_fwd` == "PERFECT" & `TYPE-BC1_rev` == "PERFECT" & n_common == 1 ~ intersect_bc,
          # Missing one barcode (fwd or rev) - and the only one barcode with the lowest lv dist
          is.na(`TYPE-BC1_fwd`) & `TYPE-BC1_rev` %in% c("PERFECT", "MISMATCHES") & `n_same_dist-BC1_rev` == 1 ~ `best_index-BC1_rev`,
          `TYPE-BC1_fwd` %in% c("PERFECT", "MISMATCHES") & is.na(`TYPE-BC1_rev`) & `n_same_dist-BC1_fwd` == 1 ~ `best_index-BC1_fwd`,
          # Different bc both ends
          `TYPE-BC1_fwd` == "PERFECT" & `TYPE-BC1_rev` == "PERFECT" & n_common != 1 ~ as.list("Chimeric"),
          # Mismatches in barcodes but they match
          `TYPE-BC1_fwd` == "MISMATCHES" & `TYPE-BC1_rev` == "MISMATCHES" & n_common == 1 ~ intersect_bc,
          `TYPE-BC1_fwd` == "PERFECT" & `TYPE-BC1_rev` == "MISMATCHES" & n_common == 1 ~ intersect_bc,
          `TYPE-BC1_fwd` == "MISMATCHES" & `TYPE-BC1_rev` == "PERFECT" & n_common == 1 ~ intersect_bc,
          # The barcodes do not match but one is perfect and the other one contains mistakes - default to the perfect one
          `TYPE-BC1_fwd` == "PERFECT" & `TYPE-BC1_rev` == "MISMATCHES" & n_common != 1 ~ `best_index-BC1_fwd`,
          `TYPE-BC1_fwd` == "MISMATCHES" & `TYPE-BC1_rev` == "PERFECT" & n_common != 1 ~ `best_index-BC1_rev`,
          # No barcode found with the right lv dist
          TRUE ~ as.list("Undetermined")
        )
    ) %>% unnest(cols = c(best_index), keep_empty = TRUE) %>%
    dplyr::rename(
      "BCtype_fwd" = `TYPE-BC1_fwd`,
      "BCtype_rev" = `TYPE-BC1_rev`,
      "bestMatches_fwd" = `best_index-BC1_fwd`,
      "bestIndexes_rev" =  `best_index-BC1_rev`,
      "NbestMatches_fwd" = `n_same_dist-BC1_fwd`,
      "NbestMatches_rev" = `n_same_dist-BC1_rev`,
      "Intersect_RevFwd" = intersect_bc,
      "NBCcommon_RevFwd" =  n_common,
      "index_unique" = best_index
    ) %>%
    mutate(EndCombination =  ifelse( is.na(BCtype_rev) | is.na(BCtype_fwd) , "BC", "BC-BC"),
           EndCombination = ifelse(index_unique == "Undetermined", NA, EndCombination),
           EndType = str_replace(paste0(BCtype_fwd, "_", BCtype_rev), "_NA|NA_", ""))


  # 7. Stats
  if(!silent){
    bc_counts = nrow(filter(real_bc, !index_unique %in% c("Undetermined", "Chimeric")))
    message(paste0("Detected Indexes in Whitelist: ", bc_counts, " (", round(100*bc_counts/barcoded_reads,2), "%)"))

    one_end = nrow(filter(real_bc, EndCombination == "BC" & !index_unique %in% c("Undetermined", "Chimeric")))
    message(paste0("Index #1 at one end: ", one_end, " (", round(100*one_end/total_reads,2), "%)"))

    two_ends = nrow(filter(real_bc, EndCombination == "BC-BC" & !index_unique %in% c("Undetermined", "Chimeric")))
    message(paste0("Index #1 at both ends: ", two_ends, " (", round(100*two_ends/total_reads,2), "%)"))

    two_ends_chimeric = nrow(filter(real_bc, index_unique == "Chimeric"))
    message(paste0("Different Index #1 between 5'/3' (Chimeric left?): ", two_ends_chimeric, " (", round(100*two_ends_chimeric/total_reads,2), "%)"))

    undetermined_stats = nrow(filter(real_bc, index_unique == "Undetermined"))
    message(paste0("Undetermined: ", undetermined_stats, " (", round(100*undetermined_stats/total_reads,2), "%)"))

  }

  # 8. Exclude barcodes supported by less than MIN_READS
  if(!silent){message(paste0("Exclude barcodes supported by less than ", MIN_READS, " reads."))}

  index <- table(real_bc$index_unique) > MIN_READS
  real_bc <- dplyr::filter(real_bc, index_unique %in% names(index)[index == TRUE])

  # 9. Regroup results with the vsearch results.
  write.table(select(real_bc, read_id, index_unique, EndCombination, EndType), paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_detected_barcodes.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(bc_sequences, paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_detected_barcodes.unfiltered.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

  # 10. Prepare the report.
  write.table(file = paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_extractBarcodes.report.txt"), sep="\t", row.names = FALSE, quote = FALSE,
              bind_cols(
                data.frame(
                  ID = OUTPUT_SUFFIX,
                  n_reads_with_vsearch_hit = barcoded_reads,
                  n_reads_correct_bc_total = nrow(filter(real_bc, !index_unique %in% c("Undetermined", "Chimeric"))),
                  n_reads_bc_oneEnd_only = nrow(filter(real_bc, EndCombination == "BC" & !index_unique %in% c("Undetermined", "Chimeric"))),
                  n_reads_bc_same_twoEnds = nrow(filter(real_bc, EndCombination == "BC-BC" & !index_unique %in% c("Undetermined", "Chimeric"))),
                  n_reads_bc_diff_twoEnds = nrow(filter(real_bc, index_unique == "Chimeric")),
                  n_reads_correct_perfect_match = sum(real_bc$EndType %in% c("PERFECT", "PERFECT_PERFECT")),
                  n_reads_correct_partial_match = sum(real_bc$EndType %in% c("MISMATCHES", "MISMATCHES_PERFECT", "PERFECT_MISMATCHES"))),
                bc_sequences_tocheck %>%
                  unnest(c(min_lv, n_same_dist)) %>%
                  group_by(min_lv, n_same_dist) %>%
                  summarise(n = n()) %>%
                  ungroup() %>%
                  dplyr::rename("min_lv_dist_bc" = "min_lv", "N_bc_same_lv_dist" = "n_same_dist"))
  )

  return("Finished")

}
