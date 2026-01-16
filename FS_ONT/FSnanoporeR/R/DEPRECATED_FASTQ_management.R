

#' Write Chimeric Segments
#'
#' REQUIRE: seqkit
#'
#' Take a FASTQ fle and divide into non_chimeric reads and chimeric ones. Chimeric ones are further spitted into segments.
#' The read_id will be modified, keeping only the prefix and appending:
#'
#' .chimeric_reads.fq.gz ==> prefix + ""
#'
#' .non_chimeric_reads.fq.gz ==> prefix + "_chimeric=FALSE_1"
#'
#' .chimeric_segments.fq.gz ==> prefix + "_chimeric=TRUE_" + read_segment_id (unique)
#'
#' @param FASTQ character. "/path/to/the/input/fastq.gz". Gziped (.gz) or not.
#' @param OUTPUT_PATH character. "/path/to/the/output/folder/"
#' @param OUTPUT_SUFFIX character. Sample ID
#' @param CHIMERIC_SEGMENTS_PATH character. Path to the blast search object created by 'detect_and_split_chimeric()'
#' @param CHIMERIC_READ_IDs character. Path to the chimeric read ids.
#'
#' @import tidyverse
#' @import ShortRead
#' @import Biostrings
#'
#' @returns
#' Three fastq files
#'
#' .chimeric_reads.fq.gz: Chimeric reads before splitting
#'
#' .non_chimeric_reads.fq.gz: Non Chimeric reads
#'
#' .chimeric_segments.fq.gz: Chimeric reads split into segments
#'
#' @export
split_and_write_fastq <- function(FASTQ = NULL,
                                  CHIMERIC_SEGMENTS_PATH = NULL,
                                  OUTPUT_PATH = NULL,
                                  CHIMERIC_READ_IDs = NULL,
                                  OUTPUT_SUFFIX = NULL){


  # 1. Split FASTQ (chimeric vs non-chimeric)
  system(intern = FALSE, ignore.stderr = TRUE,
         command = paste0(
           "seqkit grep -f ", CHIMERIC_READ_IDs, " ", FASTQ, " -o ", OUTPUT_PATH, "/", OUTPUT_SUFFIX, ".chimeric_reads.fq.gz"
           )
  )
  system(intern = FALSE, ignore.stderr = TRUE,
         command = paste0(
           "seqkit grep -vf ", CHIMERIC_READ_IDs, " ", FASTQ, " -o ", OUTPUT_PATH, "/", OUTPUT_SUFFIX, ".non_chimeric_reads.fq.gz"
           )
  )

  # 2. Open chimeric FASTQ
  fq.opened <- ShortRead::readFastq(paste0(OUTPUT_PATH, "/", OUTPUT_SUFFIX, ".chimeric_reads.fq.gz"))
  fq.ids <- stringr::str_extract(as.character(fq.opened@id), pattern = "[^ ]+")

  # 2. Match FASTQ names to the PCR_ADAPTER
  segments <- read_tsv(CHIMERIC_SEGMENTS_PATH, show_col_types = FALSE)
  segments <- filter(segments, putative_chimeric == TRUE)

  segments$fastq_id <- as.character(fq.opened@id[match(segments$read_id, fq.ids)])
  segments$read_seq <- as.character(fq.opened@sread[match(segments$read_id, fq.ids)])
  segments$read_qual <- as.character(fq.opened@quality@quality[match(segments$read_id, fq.ids)])

  # 3. Split chimeric reads into segments
  suppressWarnings(invisible(file.remove(paste0(paste0(OUTPUT_PATH, "/", OUTPUT_SUFFIX, ".chimeric_segments.fq.gz")))))

  if(any(segments$putative_chimeric)){

    ShortRead::writeFastq(mode = "w", file = paste0(OUTPUT_PATH, "/", OUTPUT_SUFFIX, ".chimeric_segments.fq.gz"),
                          ShortRead::ShortReadQ(
                            sread = Biostrings::DNAStringSet(str_sub(segments$read_seq, start = segments$segment_start, end = segments$segment_end)),
                            quality = Biostrings::BStringSet(str_sub(segments$read_qual, start = segments$segment_start, end = segments$segment_end)),
                            id = Biostrings::BStringSet(segments$fastq_id_new)
                          )
    )


  }

  return("Finished")

}



#' @title Get FASTQ stats
#'
#' @description
#' Get read_id, read_length, read_average_bp_qual (phred). Stream over a fastq file by 10000 reads to avoid overloading the memory
#'
#' @param fq character. FASTQ path
#'
#' @export
stream_FASTQ_stats <- function(fq = NULL){

  fq.opened <- ShortRead::FastqStreamer(fq, 4000)
  df <- data.frame()
  while (length(fq.subset <- ShortRead::yield(fq.opened))) {
    df <- rbind(
      data.frame(
        read_id = str_split(as.character(ShortRead::id(fq.subset)), pattern = " ", simplify = TRUE)[,1],
        read_length = ShortRead::width(fq.subset),
        average_bp_score = ShortRead::alphabetScore(Biostrings::quality(fq.subset))/ShortRead::width(fq.subset)
      ), df
    )
  }
  close(fq.opened)

  return(df)
}



