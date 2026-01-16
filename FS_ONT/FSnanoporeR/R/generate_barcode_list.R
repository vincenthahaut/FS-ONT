#' generate_index_seq_V2
#'
#' @description
#' Create a set of random barcodes.
#'
#' @param iterations numeric. Barcode length (max 14).
#' @param n numeric. Maximum number of substitutions / distance. See DNABarcodes::create.dnabarcodes.
#' @param l character. "seqlev" or "hamming"
#' @param n_bc numeric. How many BC to pick?
#' @param nucl character. Tried 'conway' (fast) and ashlock (slow).
#'
#' @importFrom reshape2 melt
#' @import tidyverse
#' @import stringdist
#'
#' @return
#' data.frame with the following columns:
#' barcodes sequence
#' minimal levenstein distance of this barcode to the other in that set
#' set (=iteration)
#' min.lv.set minimal levenstein distance of that set
#'
#' @examples
#' bc18 <- randomBC(iterations = 500, n = 200000, l = 18, nucl = c("A", "T", "C", "G"))
#' bc19 <- randomBC(iterations = 500, n = 200000, l = 19, nucl = c("A", "T", "C", "G"))
#' bc20 <- randomBC(iterations = 500, n = 200000, l = 20, nucl = c("A", "T", "C", "G"))
#' bc21 <- randomBC(iterations = 500, n = 200000, l = 21, nucl = c("A", "T", "C", "G"))
#' bc22 <- randomBC(iterations = 500, n = 200000, l = 22, nucl = c("A", "T", "C", "G"))
#' # Select bc20 with the least number of BC harboring a minimal lv distance to the other barcode (=6)
#' index <- bc20 %>%
#'   group_by(it_n) %>%
#'   summarise(bc_lv6 = sum(min_dist == 6, na.rm = T))
#'
#' final.bc <- bc20 %>%
#'   filter(it_n == index$it_n[min(index$bc_lv6)]) %>%
#'   filter(min_dist >= 7) %>%
#'   head(384)
#'
#' hist(lv_dist_to_df(final.bc$row, final.bc$row)$min_dist)
#' hist(str_count(final.bc$row, "G|C")/20)
#'
#' final.bc <- final.bc$row
#' names(final.bc) <- paste0("BC", 1:384)
#'
#' ShortRead::writeFasta(Biostrings::DNAStringSet(final.bc), "/home/vincent.hahaut/Desktop/IOB/6_Nanopore-FS/FSnanoporeR/inst/extdata/whitelist.bc20.fa")
#'
#' @export
randomBC <- function(iterations = 100, n = 100000, l = 20, n_bc = 96, nucl = c("A", "T", "C", "G")){

  library(tidyverse)

  lv_dist_to_df <- function(barcodes_1 = NULL, barcodes_2 = NULL, method = "lv"){
    lv_dist <- stringdist::stringdistmatrix(barcodes_1, barcodes_2, method = method) %>% as.data.frame()
    colnames(lv_dist) <- barcodes_2
    row.names(lv_dist) <- barcodes_1
    lower.tri <- lower.tri(lv_dist) %>% as.data.frame() %>%
      rownames_to_column("row") %>%
      gather(., key = "col", value = "value", -row)
    lv_dist <- lv_dist %>%
      rownames_to_column("row") %>%
      gather(., key = "col", value = "value", -row)


    reshape2::melt(lv_dist, varnames = c("row", "col"))

    lv_dist %>%
      dplyr::filter(lower.tri$value) %>%
      group_by(row) %>%
      summarise(min_dist = min(value)) %>%
      ungroup() %>%
      # Put back the missing element (first)
      bind_rows(data_frame(row = barcodes_1[1], min_dist = NA)) %>%
      return()
  }

  set.seed(1000)
  df.l <- data_frame()

  # 1. Generate X random sequences of length l
  seqs <- unique(sapply(1:n, function(x) paste0(sample(nucl, l, replace = TRUE), collapse = "")))

  # 2. Remove BC with trimer or starting/ending with the same base as the adjacent one in the barcoding primer
  # Remove barcodes with outlier GC content
  seqs <- seqs[!str_detect(seqs, "AAA|TTT|CCC|GGG|A$|^TT")]
  seqs <- seqs[str_count(seqs, "G|C")/l > 0.35 & str_count(seqs, "G|C")/l < 0.65]

  # 3. For i interations, select X random ones and calculate their minimal levenstein distance
  for(i in 1:iterations){
    seqs.tmp <- sample(seqs, n_bc, replace = FALSE)
    df.l <- bind_rows(df.l,
                      lv_dist_to_df(seqs.tmp, seqs.tmp, method = "lv") %>%
                        mutate(bc_length = l,
                               it_n = i,
                               min.lv.set = min(min_dist, na.rm= T))
    )
  }
  return(df.l)

}


########

# FIRST ITERATION OF THE BARCODE SET

#' generate_index_seq_V1
#'
#' @description
#' Create a set of random barcodes using DNABarcodes. Can be very slow if long barcodes are required (>24hrs).
#' This function is tailored for this particular application and was not extensively benchmarked with other parameters than those proposed.
#'
#' @param size numeric. Barcode length (max 14).
#' @param dist numeric. Maximum number of substitutions / distance. See DNABarcodes::create.dnabarcodes.
#' @param metric character. "seqlev" or "hamming"
#' @param heuristic character. Tried 'conway' (fast) and ashlock (slow).
#' @param outpath character. /path/to/results/
#' @param threads numeric. DNABarcodes::create.dnabarcodes threads. Only useful with heuristic = ashlock.
#' @param context character. DNA barcode final context marked with Ns "CGACGCTCTTCCGATCTNNNNNNNNNNNNNAAGCAGTGGTATCAACGCA".
#' @param nanopore_bc character. "/path/to/the/output/native/barcoding/barcode/sequence.fasta
#' @param nanopore_kit_prefix character. Example: "NB"
#'
#' @importFrom DNABarcodes create.dnabarcodes
#' @import stringdist
#' @import ShortRead
#'
#' @details
#' 1. Run DNABarcodes::create.dnabarcodes using the defined parameters.
#' ==> Idealy I aimed at 12bp bc / seqlev=4 but this function kept returning 320 bc (60 short from full-plate).
#' ==> Went for 13bp instead.
#' 2. Filter out sequences with triplets (including those with the barcode edges, ex. T .. A)
#' 3. Filter out sequences with >2 homodimers (including those at the barcode edge, ex. T .. A)
#' 4. Filter out GC content >65% and <35%.
#' 5. Exclude barcodes too close from the nanopore_bc / adapters (lv dist <10).
#' 6. Sample randomly 384 barcodes 10000 times and calculate their levenstein distance.
#' 7. Get the barcodes list with the smallest mean levenstein distance to maximise the chances during demultiplexing.
#' 8. Return the barcodes as a fasta file.
#'
#'@export

# Create Barcodes
# Can take many hours to run!
generate_index_seq <- function(size = 13,
                               dist = 4,
                               metric = "seqlev",
                               threads = 24,
                               heuristic = "conway",
                               outpath = "Desktop/",
                               context = "CGACGCTCTTCCGATCTNNNNNNNNNNNNNAAGCAGTGGTATCAACGCA",
                               nanopore_bc = "/opt/ont/guppy/data/barcoding/barcodes_masked.fasta",
                               nanopore_kit_prefix = "NB"){

  #suppressPackageStartupMessages(library(DNABarcodes))
  #suppressPackageStartupMessages(library(stringdist))
  #suppressPackageStartupMessages(library(ShortRead))
  #suppressPackageStartupMessages(library(tidyverse))

  out <- paste0(outpath, "/Barcodes_seq", size, "_seql", dist, "_heuristic", heuristic, "_metric", metric, ".txt")
  if(!file.exists(out)){
    message("Takes >24hrs")
    bc <- DNABarcodes::create.dnabarcodes(size,
                                          dist = dist,
                                          metric=metric,
                                          cores=threads,
                                          heuristic = heuristic,
                                          filter.triplets = TRUE,
                                          filter.gc = FALSE,
                                          population= ifelse(heuristic == "ashlock", 500, 1),
                                          filter.self_complementary = FALSE,
                                          # didn't see any benefit in going above
                                          iterations=1)
    write.table(bc, out, quote = FALSE, row.names = FALSE, col.names = FALSE)
    message(paste0("Creadted: ", length(bc), " Barcodes"))
  } else {
    bc <- read.table(out, sep = "\t", header = FALSE)
    message(paste0("Loaded: ", nrow(bc), " Barcodes"))
  }

  # 2.1. Flag Homopolymers
  homodimer <- sapply(bc, function(x) str_count(x, pattern = "AA|TT|GG|CC"))
  homotrimer <- sapply(bc, function(x) str_detect(x, pattern = "AAA|TTT|GGG|CCC"))

  # 2.2. Flag context creating more homopolymers
  bc_pos <- str_locate(context, "N.*N")
  downstream <- str_sub(context, bc_pos[,1]-1, bc_pos[,1]-1)
  upstream <- str_sub(context, bc_pos[,2]+1, bc_pos[,2]+1)

  downstream_homodimer <- sapply(bc$V1, function(x) str_detect(x, pattern = paste0("^", downstream)))
  upstream_homodimer <- sapply(bc$V1, function(x) str_detect(x, pattern = paste0(upstream, "$")))

  downstream_homotrimer <- sapply(bc$V1, function(x) str_detect(x, pattern = paste0("^", paste0(rep(downstream,2), collapse = ""))))
  upstream_homotrimer <- sapply(bc$V1, function(x) str_detect(x, pattern = paste0(paste0(rep(upstream,2),  collapse = ""),"$")))


  # 2.3. Summarise & filter
  # Max 3 homodimers, including the edge nucleotides.
  bc <- data_frame(
    bc = bc$V1,
    hd = homodimer[,1],
    ht = homotrimer[,1],
    down_hd = downstream_homodimer,
    up_hd = upstream_homodimer,
    down_htri = downstream_homotrimer,
    up_htri = upstream_homotrimer
  ) %>%
    mutate(GC = 100*(str_count(bc, "G|C"))/size,
           total_hd = hd + down_hd + up_hd) %>%
    # Specific to the ISPCR
    dplyr::filter(!str_detect(pattern = "A$", bc)) %>%
    dplyr::filter(ht != TRUE & down_htri != TRUE & up_htri != TRUE) %>%
    dplyr::filter(total_hd < 3) %>%
    dplyr::filter(GC > 35 & GC < 65)

  # 3. Check distances
  # Use levenstein as it is faster than seqlev and provide similar results here

  # 3.1. Verify against nanopore sequences
  lv_dist_to_df <- function(barcodes_1 = NULL, barcodes_2 = NULL, method = "lv"){
    lv_dist <- stringdist::stringdistmatrix(barcodes_1, barcodes_2, method = method) %>% as.data.frame()
    colnames(lv_dist) <- barcodes_2
    row.names(lv_dist) <- barcodes_1
    lower.tri <- lower.tri(lv_dist) %>% as.data.frame() %>%
      rownames_to_column("row") %>%
      gather(., key = "col", value = "value", -row)
    lv_dist <- lv_dist %>%
      rownames_to_column("row") %>%
      gather(., key = "col", value = "value", -row)


      reshape2::melt(lv_dist, varnames = c("row", "col"))

    lv_dist %>%
      dplyr::filter(lower.tri$value) %>%
      group_by(row) %>%
      summarise(min_dist = min(value)) %>%
      ungroup() %>%
      # Put back the missing element (first)
      bind_rows(data_frame(row = barcodes_1[1], min_dist = NA)) %>%
      return()
  }

  ont_bc <- ShortRead::readFasta(nanopore_bc)
  ont_bc <- data.frame(bc = ont_bc@sread, bc_id = ont_bc@id) %>%
    dplyr::filter(grepl(paste0("^", nanopore_kit_prefix), bc_id))

  lv_dist_ont <- lv_dist_to_df(barcodes_1 = bc$bc, barcodes_2 = ont_bc$bc)

  bc <- mutate(bc, min_lv_ont_bc = lv_dist_ont$min_dist) %>%
    dplyr::filter(min_lv_ont_bc > 10)

  # 4.2. Maximise the bc-bc distance
  # Randomly sample 384 barcodes and calculate their distance.
  # Select the set with the lowest average.
  set.seed(42)
  bc.sampling.l <- list()
  bc.sampling.df <- data.frame()
  attempts <- 10000
  n_bc = 384

  bc.optimisation.results <- lapply(1:attempts, function(i) {
    message(i)
    bc.sampling <- sample_n(bc, n_bc, replace = FALSE)
    lv_dist_to_df(bc.sampling$bc, bc.sampling$bc) %>%
      mutate(
        attempt = i,
        mean_dist = mean(min_dist, na.rm = T),
        dist_exact = sum(min_dist == dist, na.rm = T),
        dist_p1 = sum(min_dist == dist+1, na.rm = T),
        dist_p2 = sum(min_dist == dist+2, na.rm = T),
        dist_p3 = sum(min_dist == dist+3, na.rm = T),
      )
  })

  # 6. Extract the best barcode set
  best_set <- sapply(bc.optimisation.results, function(x) unique(x["mean_dist"]))
  best_set <- unlist(best_set)
  best_set <- which(best_set == min(best_set))
  bc.final <- right_join(bc, bc.optimisation.results[[best_set]], by = c("bc" = "row"))

  # 7. Save the results + info
  out <- paste0(outpath, "/Barcodes_seq", size, "_seql", dist, "_heuristic", heuristic, "_metric", metric, "_selectedSet.txt")
  write.table(bc.final, out, quote = FALSE, row.names = FALSE, col.names = FALSE)

  # 8. Save the results as a fasta file
  barcodes.final.fa <- bc.final$bc
  names(barcodes.final.fa) <- paste0("BCFS_", 1:n_bc)

  out <- paste0(outpath, "/ONT_Barcodes_FS.fa")
  barcodes.final <- DNAStringSet(barcodes.final.fa, use.names = TRUE)
  ShortRead::writeFasta(barcodes.final, file = out)

  return(bc.final)

}


