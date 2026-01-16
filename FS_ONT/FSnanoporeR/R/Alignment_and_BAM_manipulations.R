#' run_minimap2
#'
#' @description
#' Map FASTQ file using Minimap2 and filter out unmapped reads.
#'
#' Requires minimap2 path. Requires samtools installed.
#'
#' @param FASTQ character. "/path/to/the/input/fastq.gz". Gziped (expect '.gz') or not. Accept also a vector of FASTQ file paths.
#' @param OUTPUT_SUFFIX character. Sample prefix id.
#' @param MINIMAP2_Binaries character. "/path/to/minimap2" binaries
#' @param MINIMAP2_REF character. "/path/to/minimap2_reference.fa" FASTA or .mni. Using FASTA will significantly slow down the process.
#' @param MINIMAP2_BED character. "/path/to/minimap2_gene_reference.bed" prepared using paftools.js
#' @param OUTPUT_DIR character. "/path/to/the/output/".
#' @param ISOQUANT_MODE boolean. If TRUE, use the -Y flag to mark secondary alignments as soft-cliping.
#' @param CHIMERIC_ON boolean. If TRUE will select sam flag '260' else '2308'.
#' @param STRAND_SPECIFIC boolean. Should it be run only considering the forward strand (-uf).
#' @param THREADS numeric. Number of threads.
#'
#' @details
#' https://github.com/lh3/minimap2/blob/master/misc/README.md
#'
#' @returns Creates a sorted BAM file at paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_sorted.bam").
#'
#' @export
run_minimap2 <- function(FASTQ = NULL,
                         OUTPUT_SUFFIX = NULL,
                         MINIMAP2_Binaries = NULL,
                         MINIMAP2_REF = NULL,
                         MINIMAP2_BED = NULL,
                         OUTPUT_DIR = "./",
                         CHIMERIC_ON = FALSE,
                         STRAND_SPECIFIC = FALSE,
                         ISOQUANT_MODE = FALSE,
                         THREADS = 1){

  if(!file.exists(MINIMAP2_Binaries) & MINIMAP2_Binaries != "minimap2"){
    stop("Issues with minimap2-binary")
  }
  if(!file.exists(MINIMAP2_REF) | !file.exists(MINIMAP2_BED)){
    stop("Issues with minimap2-related paths")
  }

  if(grepl("not found", system(intern = TRUE, "type samtools &>/dev/null"))){
    stop("Please install samtools")
  }

  if(CHIMERIC_ON){
    sam_FLAG <- 260
  } else {
    sam_FLAG <- 2308
  }

  if(length(FASTQ) > 1){
    FASTQ <- paste0(FASTQ, collapse = " ")
  }

  # 1. Mapping
  system(intern = FALSE,
         command = paste0(
           MINIMAP2_Binaries, " ",
           ifelse(STRAND_SPECIFIC == TRUE, "-uf ", "-ub "),
           ifelse(ISOQUANT_MODE == TRUE, "-Y ", " "),
           "-ax splice --MD ",
           "-t ", THREADS, " ",
           "--junc-bed ", MINIMAP2_BED, " ",
           MINIMAP2_REF, " ",
           FASTQ,
           "| samtools view -F ", sam_FLAG, " -b",
           "| samtools sort -@ ", THREADS, " -o", paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_sorted.bam"), " >/dev/null 2>&1")
  )

  # 3. Indexing
  system(intern = FALSE,
         command = paste(
           "samtools index",
           paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_sorted.bam"))
  )

  return("Finished")

}

#' load_bam
#'
#' @description
#' Load a BAM file to R.
#' If getSoftClipping = TRUE, will also report the position on the read of left/right-clippings and mapping.
#'
#' @param BAM_path character. "/path/to/the/input/file.sorted.bam". Requires an indexed file in the same folder.
#' @param getSoftClipping boolean. If TRUE, will report the right / left clipping (soft & hard) start/end/width in the read. Useful to spot putative chimeric reads and compare read mapping segments.
#' @param clip_flag_width numeric. Flag as "putative chimeric" reads with a soft-clipped / hard-clipped value superior to 'soft_clip_flag_width' (left/right).
#' @param FIELDS character vector. BAM fields to read. Expect the following for chimeric detection: c("qname","flag", "rname","strand","pos", "qwidth", "mapq", "cigar")
#'
#' @import tidyverse
#' @import Rsamtools
#' @import S4Vectors
#' @importFrom GenomicAlignments cigarWidthAlongReferenceSpace
#' @importFrom GenomicAlignments cigarRangesAlongQuerySpace
#' @importFrom parallel mclapply
#'
#' @returns BAM file to data.frame.
#'
#' @export
load_bam <- function(BAM_path = NULL, tidy_chimeric = TRUE, clip_flag_width = 300, THREADS = 1, FIELDS = c("qname","flag", "rname","strand","pos", "qwidth", "mapq", "cigar")){

  #suppressPackageStartupMessages(library(Rsamtools))
  #suppressPackageStartupMessages(library(tidyverse))
  #suppressPackageStartupMessages(library(GenomicAlignments))

  params <- Rsamtools::ScanBamParam(what=FIELDS)
  bam <- Rsamtools::scanBam(BAM_path, param = params)

  if(length(bam[[1]]$qname) == 0){

    message("No read detected!")

  } else {

    mapping_width <- GenomicAlignments::cigarWidthAlongReferenceSpace(bam[[1]]$cigar)

    bam <- data.frame(
      read_id = bam[[1]]$qname,
      sam_flag = bam[[1]]$flag,
      seqnames = bam[[1]]$rname,
      map_start = bam[[1]]$pos,
      map_end = bam[[1]]$pos + mapping_width-1,
      qwidth = bam[[1]]$qwidth,
      mapq = bam[[1]]$mapq,
      strand = bam[[1]]$strand,
      cigar = bam[[1]]$cigar) %>%
      group_by(read_id) %>%
      # Unique ID to account for supplementary alignments
      mutate(mapped_read_id = paste0(read_id, ".", 1:n()),
             genome_position = paste0(seqnames, ":", map_start, "-", map_end, ":", strand)) %>%
      ungroup()

    if(tidy_chimeric == TRUE){

      # 0. Update read_ids
      # Looking for chimeric reads rely on the defining the chromosome as read_id+seqnames
      # Unfortunately, seqnames can become too long for GRanges (overflow error)
      # Convert read IDs to a simpler format
      bam$simplified_read_id <- paste0("mapread_", as.numeric(factor(str_replace(bam$read_id, "_chim.*$", ""))))
      bam$simplified_read_id <- ifelse(grepl("_chim", bam$read_id),
                                       paste0(bam$simplified_read_id, "_chim_", str_extract(bam$read_id, "[^_]+$")),
                                       bam$simplified_read_id)

      # 1. Extract relevant information from the BAM file

      # 1.1. Read Length
      read_length <- parallel::mclapply(bam$cigar, mc.cores = THREADS, function(x)
        GenomicAlignments::cigarRangesAlongQuerySpace(x, before.hard.clipping = TRUE)[[1]] %>% S4Vectors::end(.) %>% max) %>%
        unlist

      # 1.2. Read info
      bam.info <- select(bam, mapped_read_id, simplified_read_id, read_id, strand, sam_flag, genome_position, mapq, qwidth, seqnames, map_start, map_end) %>%
        mutate(sam_flag = ifelse(sam_flag %in% c(0,16), "PRIMARY", "SUPPLEMENTARY"),
               read_length = read_length) %>%
        dplyr::rename("read_orientation" = strand,
                      "mapping_width" = qwidth)

      # 1.3. Soft-/hard-clipping start-end in the read
      clipping <- GenomicAlignments::cigarRangesAlongQuerySpace(bam$cigar, ops=c("S", "H"), before.hard.clipping = TRUE, with.ops = TRUE)
      names(clipping) <- bam$mapped_read_id
      clipping.df <- as.data.frame(clipping) %>%
        select(group_name, names, start, end, width)
      colnames(clipping.df) <- c("mapped_read_id", "clip_type" , "read_start", "read_end", "width")

      # 1.4. Mapped portion of the read (start-end)
      map_pos <- GenomicAlignments::cigarRangesAlongQuerySpace(bam$cigar, with.ops = TRUE, ops = c("M","N","I","D"), before.hard.clipping = TRUE)
      names(map_pos) <- bam$mapped_read_id
      map_pos <- parallel::mclapply(map_pos, mc.cores = THREADS, function(x) data.frame(read_start = min(S4Vectors::start(x)), read_end = max(S4Vectors::end(x))))
      map_pos.df <- do.call(rbind, map_pos) %>%
        rownames_to_column("mapped_read_id") %>%
        mutate(clip_type = "M",
               width = read_end - read_start)

      # 2. Combine BAM / Read-position data
      cigar.info <- bind_rows(map_pos.df, clipping.df) %>%
        # Determine if the read portion is mapped or clipped (left or right-end of the read)
        mutate(clip_pos =
                 case_when(
                   clip_type == "M" ~ "Mapped",
                   read_start == 1 & clip_type != "M" ~ "left_clip",
                   TRUE ~ "right_clip"),
        # Flag putative chimeric
        # Defined as having a large clipping (left or right)
               putative_chimeric_mapping = width > clip_flag_width & clip_type != "M")

      cigar.info <- left_join(cigar.info, bam.info, by = c("mapped_read_id")) %>%
        mutate(
          read_id = str_replace(read_id, "\\..$", ""),
          read_group = str_extract(read_id, "(?<=_)\\d+$")
        )

      # 3. Correct read orientation
      # CIGAR are not reported 5'-3' but relative to the read orientation.
      cigar.info <- cigar.info %>%
        mutate(
          read_start_corrected = ifelse(read_orientation == "+", read_start, read_length-read_end),
          read_end_corrected = ifelse(read_orientation == "+", read_end, read_length-read_start)
        ) %>%
        arrange(mapped_read_id, read_start_corrected, read_end_corrected, read_group)


      if(any(grepl("_chimeric=", cigar.info$read_id))){
        cigar.info$chimera_flag = str_extract("(?<=_chimeric=).*[^_]+(?=_)", string = cigar.info$read_id)
      } else {
        cigar.info$chimera_flag = NA_character_
      }

      return(cigar.info)

    } else {
      return(bam)
    }

  }

}
