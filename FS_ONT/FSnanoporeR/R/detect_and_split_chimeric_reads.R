#' @title Run BLAST using a single target sequence
#'
#' @description
#' Convert FASTQ into FASTA. Create a blast 5 database using the adapter sequence and map the reads against it.
#'
#' @param FASTQ character. "/path/to/the/input/fastq.gz". Gziped (expect '.gz') or not.
#' @param OUTPUT_DIR character. "/path/to/output/"
#' @param OUTPUT_SUFFIX character.
#' @param PCR_ADAPTER_REV character. Sequence searched with BLASTn to define the read orientation and chimerism.
#' @param PCR_ADAPTER_FWD character. Sequence searched with BLASTn to define the read orientation and chimerism.
#' @param THREADS numeric.
#' @param IDENTITY numeric. Percentage of identity between query and "PCR_ADAPTER" to be considered as a valid "PCR_ADAPTER" hit.
#' @param silent boolean. If TRUE will print various messages.
#'
#' @return
#' data.frame - _blast_search.txt
#'
#' @export
run_blast <- function(FASTQ = NULL,
                      OUTPUT_DIR = NULL,
                      OUTPUT_SUFFIX = NULL,
                      PCR_ADAPTER_REV = "AAGCAGTGGTATCAACGCAGAGT",
                      PCR_ADAPTER_FWD = "AAGCAGTGGTATCAACGCAGAGT",
                      THREADS = 1,
                      IDENTITY = 65,
                      silent = FALSE){

  silent = as.logical(toupper(silent))

  if(grepl("not found", system(intern = TRUE, "type blastn &>/dev/null"))){
    stop("Please install blast")
  }

  # 1. Convert FASTQ to FASTA
  if(!silent){message("run_blast: Convert FASTQ to FASTA")}

  system(intern = FALSE, ignore.stdout = silent, ignore.stderr = silent, paste0(
    "seqkit fq2fa --threads ", THREADS, " ", paste0(FASTQ, collapse = " "), " -o ", OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_tmp.fa"))

  # 2. Create BLAST search database
  if(!silent){message("run_blast: Create the BLAST database")}
  blast_db <- paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_ADAPTER_rblast_db.tmp.fa")
  if(file.exists(blast_db)){file.remove(blast_db)}
  write(">PCR_ADAPTER_FWD", blast_db)
  write(PCR_ADAPTER_FWD, blast_db, append = TRUE)
  write(">PCR_ADAPTER_REV", blast_db, append = TRUE)
  write(as.character(ShortRead::reverseComplement(Biostrings::DNAStringSet(PCR_ADAPTER_REV))), blast_db, append = TRUE)

  system(intern = FALSE, ignore.stdout = TRUE, ignore.stderr = FALSE,
         command = paste(
           "makeblastdb -in", blast_db,
           "-parse_seqids -dbtype nucl")
  )

  # 3. BLAST Search
  # blast -perc_identity = min. identity percentage allowed.
  # Typically with nchar(PCR_ADAPTER) + these settings the minimal returned identity is slightly above 70% anyway.
  # Looking for the full ADAPTER2[N..N]ADAPTER1 does not work well with BLAST.
  # BLAST does not multithread well.
  if(!silent){message("run_blast: Detect PCR Adapters using BLASTn-Short")}
  system(intern = TRUE,
         command = paste(
           "blastn  -task blastn-short",
           "-query", paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_tmp.fa"),
           "-db", blast_db,
           "-strand plus",
           "-word_size 11",
           "-outfmt '6 qseqid sseqid pident mismatch gapopen qstart qend qlen qseq qcovs sstart send sseq sstrand length evalue bitscore'",
           "-num_threads", THREADS,
           "-gapopen", 1,
           "-gapextend", 1,
           "-window_size", 0,
           "-perc_identity", IDENTITY,
           "-out", paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_blast_search.txt")
         )
  )

  # Clean-up
  if(!silent){message("run_blast: Cleanup")}
  invisible(file.remove(paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_tmp.fa")))
  invisible(unlink(paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_ADAPTER_rblast_db.tmp.fa.*$")))

}







#' @title Detect and Split Chimeric Reads
#'
#' @description
#' Summary: Identity chimeric reads and split them into segments.
#'
#' 1. Run BLASTn-Short looking for any "PCR_ADAPTER" in the read with an identity > "IDENTITY".
#'
#' 2. Chimeric reads are defined as having at least a "PCR_ADAPTER" located > "CHIMERIC_DETECTION_WINDOW" base-pairs of the 5' or 3' end.
#'
#' 3. Ideal reads have a 5' PCR_ADAPTER(+) and a 3' PCR_ADAPTER(-). Parse the BLAST output and define split the chimeric reads into segments according to the "PCR_ADAPTER" orientations such as:
#'
#'
#' Read:
#' (5')...<==|==>...<==|==>...(3')
#'
#'
#' Segments:
#'
#'
#' (+)==>...<==(-)
#'
#'
#' (+)==>...(read_end)
#'
#'
#' (read_start)...<==(-)
#'
#'
#' 4. Filter out segments <100 bp.
#'
#' 5. Add a 50bp at the "PCR_ADAPTER" 3' end to encompass the cell barcode #1 if present. Currently hard-coded to 50bp as I except than large window values may be problematic for downstream analysis. This approach will likely extract a few bases from neighboring chimeric segments. However, minimap2 will align them as soft-clipped sequences given their short size and they are unlikely to interefere with the rest of the analysis.
#'
#' 6. Correct the read length to avoid going <0 bp or >read_length.
#'
#' 7. Save the BLAST result table.
#'
#' 8. Use 'split_and_write_fastq()' to create three FASTQ.gz files: non-chimeric reads, chimeric reads and chimeric segments (=splitted reads).
#'
#' It should be noted that any method solely based on the detection of the "PCR_ADAPTER" to split chimeric reads will miss the following chimeric case:
#'
#'
#' (+) ==>... | ... <== (-)
#'
#'
#' Where two reads are ligated in opposite orientations and lack a central "PCR_ADAPTER". Although rare, these cases are handled by 'chimeric_mapping_validation()'.
#'
#' NOTE: require blast installed. Expect BLAST version compatible with blast database version 5. BLAST is used instead of VSEARCH because it can return multiple similar hits per read.
#'
#' @param FASTQ character. "/path/to/the/input/fastq.gz". Gziped (expect '.gz') or not.
#' @param OUTPUT_DIR character. "/path/to/output/"
#' @param OUTPUT_SUFFIX character.
#' @param PCR_ADAPTER_REV character. Sequence searched with BLASTn to define the read orientation and chimerism.
#' @param PCR_ADAPTER_FWD character. Sequence searched with BLASTn to define the read orientation and chimerism.
#' @param CHIMERIC_DETECTION_WINDOW numeric. Consider a read as chimeric if it has a "PCR_ADAPTER" sequence located >'CHIMERIC_DETECTION_WINDOW' bp from the 5' or 3' edge.
#' @param IDENTITY numeric. Percentage of identity between query and "PCR_ADAPTER" to be considered as a valid "PCR_ADAPTER" hit.
#' @param silent boolean. If TRUE will print various messages.
#'
#' @import tidyverse
#' @import ShortRead
#' @import Biostrings
#'
#' @return
#' data.frame - _chimericSegments.info.txt
#'
#' FASTQ.gz #1 - _chimeric_reads.fq.gz
#'
#' FASTQ.gz #2 - _chimeric_segments.fq.gz
#'
#' FASTQ.gz #3 - _non_chimeric_reads.fq.gz
#'
#' @export
split_chimeric_reads <- function(FASTQ = NULL,
                                 OUTPUT_DIR = NULL,
                                 OUTPUT_SUFFIX = NULL,
                                 PCR_ADAPTER_REV = "",
                                 PCR_ADAPTER_FWD = "",
                                 CHIMERIC_DETECTION_WINDOW = 200,
                                 IDENTITY = 65,
                                 silent = FALSE){

  if(!silent){message("1. Run BLAST")}
  # Slowest part
  run_blast(FASTQ = FASTQ,
            OUTPUT_DIR = OUTPUT_DIR,
            OUTPUT_SUFFIX = OUTPUT_SUFFIX,
            PCR_ADAPTER_REV = PCR_ADAPTER_REV,
            PCR_ADAPTER_FWD = PCR_ADAPTER_FWD,
            IDENTITY = IDENTITY,
            silent = silent)

  # 4.1. Load BLAST results
  if(!silent){message("2. Load BLAST results")}
  blast.res <- read_tsv(file = paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_blast_search.txt"),
                        show_col_types = FALSE,
                        col_names = c("read_id","sseqid","pident","mismatch","gapopen","qstart","qend","qlen","qseq","qcovs","sstart","send","sseq","sstrand","length","evalue","bitscore")) %>%
    # Exclude too low match length
    dplyr::filter(length >= floor(0.6*nchar(PCR_ADAPTER_REV)) & sseqid == "PCR_ADAPTER_REV" |
                    length >= floor(0.6*nchar(PCR_ADAPTER_FWD)) & sseqid == "PCR_ADAPTER_FWD") %>%
    # Recode match orientation
    mutate(sstrand = ifelse(sseqid == "PCR_ADAPTER_REV", "minus", sstrand)) %>%
    # Detect chimeric
    # Flag reads with more than 2 adapters as chimeric
    # Flag reads with at least one adapter too far from the edge
    mutate(putative_chimeric = qstart > CHIMERIC_DETECTION_WINDOW & qstart < (qlen - CHIMERIC_DETECTION_WINDOW)) %>%
    group_by(read_id) %>%
    mutate(n_adapters = n(),
           putative_chimeric = ifelse(n_adapters > 2, TRUE, putative_chimeric),
           putative_chimeric = any(putative_chimeric)) %>%
    ungroup()

  if(all(!blast.res$putative_chimeric)){
    # Stop function but does not crash script if missing chimeric reads
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop("No Chimeric Read Detected. Moving to next step.")
  }

  # 5. Tidy & Annotate PCR adapter orientation
  if(!silent){message("3. Split Chimeric Reads into Segments")}

  # 5.1. NON_CHIMERIC
  blast.res.segments_non_chimeric <- filter(blast.res, putative_chimeric == FALSE)

  blast.res.segments_non_chimeric$chimera_type <- "==>...<=="
  blast.res.segments_non_chimeric$chimera_type[blast.res.segments_non_chimeric$n_adapters == 1 & blast.res.segments_non_chimeric$sstrand == "plus"] <- "==>..."
  blast.res.segments_non_chimeric$chimera_type[blast.res.segments_non_chimeric$n_adapters == 1 & blast.res.segments_non_chimeric$sstrand == "minus"] <- "...<=="

  blast.res.segments_non_chimeric <- blast.res.segments_non_chimeric %>%
    select(read_id, qlen, chimera_type) %>%
    distinct() %>%
    mutate(segment_start = 1,
           putative_chimeric = FALSE,
           segment_width = qlen,
           read_length = qlen,
           n_segments = 1,
           fastq_id_new = paste0(read_id, "_chimeric=FALSE_1") ) %>%
    dplyr::rename("segment_end" = qlen) %>%
    dplyr::filter(segment_width > 100) %>%
    select(read_id, segment_start, segment_end, segment_width, read_length, chimera_type, putative_chimeric, n_segments, fastq_id_new)

  # 5.2. CHIMERIC
  win = 30

  blast.res.segments_chimeric <- filter(blast.res, putative_chimeric == TRUE) %>%
    select(read_id, qstart, qend, sstrand, qlen, putative_chimeric) %>%
    group_by(read_id) %>%
    arrange(qstart) %>%
    # Get read edges and adjust
    # Arbitrary adds 'win' bp as this would cover any barcode
    # May enter in the next read and create short duplicates of the edge sequences but should be limited to <win which will not interfere with mapping or barcode assignment.
    # Undetectable chimeric case with this method:
    # ==>(plus)... | ...(minus)<== where two reads without a detected middle ISPCR are found (false negative or absent).
    reframe(
      segment_start = dplyr::case_when(
        # (|==>...) and (|==>...<==)
        sstrand == "plus" ~ qstart-win,
        # (|==>...<==) - handle the (<==)
        sstrand == "minus" & dplyr::lag(sstrand) == "plus" ~ dplyr::lag(qstart)-win,
        # Account for (read_start...<==) and (<==|...<==)
        sstrand == "minus" & dplyr::lag(sstrand, default = "minus") != "plus" ~ dplyr::lag(qend, default = 1),
      ),
      segment_end = dplyr::case_when(
        # (...<==|) and (==>...<==|)
        sstrand == "minus" ~ qend+win,
        # (==>...<==|) - handle the (==>)
        sstrand == "plus" & lead(sstrand) == "minus" ~ lead(qend)+win,
        # Account for (==>...read_end) and (==>...|==>)
        sstrand == "plus" & lead(sstrand, default = "plus") != "minus" ~ lead(qstart, default = unique(qlen)),
      ),
      chimera_type = dplyr::case_when(
        sstrand == "plus" & lead(sstrand) == "minus" ~ "==>...<==",
        sstrand == "minus" & dplyr::lag(sstrand) == "plus" ~ "==>...<==",
        sstrand == "plus" ~ "==>...",
        sstrand == "minus" ~ "...<=="),
      segment_width = segment_end - segment_start,
      read_length = unique(qlen),
      putative_chimeric = putative_chimeric
    ) %>%
    # Collapse (plus)==>...<==(minus) cases.
    # The annotation method above double them.
    distinct() %>%
    # Filter out too short segments
    # Could happen for instance when two non-barcoded molecules are ligated
    # <== ... (<win-bp) ... ==>
    # Or if TSO invade each others (should be prevented by the 5'biotin). Found one case.
    dplyr::filter(segment_width > 100) %>%
    # Correct read positions.
    mutate(segment_start = if_else(segment_start < 1, 1, segment_start),
           segment_end = if_else(segment_end > read_length, read_length, segment_end),
           # Number of segment detected
           n_segments = n(),
           # Updated FASTQ ID
           fastq_id_new = paste0(read_id, "_chimeric=TRUE_", 1:n()) ) %>%
    ungroup()

  # 5.3. Combine
  blast.res.segments <- bind_rows(
    blast.res.segments_non_chimeric,
    blast.res.segments_chimeric)

  if(!silent){message("4. Save Chimeric Reads Infos")}
  write.table(blast.res.segments,
              file = paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_chimericSegments.info.txt"),
              quote = FALSE,
              sep = "\t",
              row.names = FALSE)

  write.table(blast.res.segments_chimeric$read_id,
              file = paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_chimericReads.ids.txt"),
              quote = FALSE,
              col.names = FALSE,
              sep = "\t",
              row.names = FALSE)

  # 5.4. Cleanup for RAM
  blast.res.segments_non_chimeric <- blast.res.segments_chimeric <- NULL

  # 6. Statistics
  if(!silent){
    n_reads = length(unique(blast.res.segments$read_id))
    n_non_chimeric = length(unique(blast.res.segments$read_id[blast.res.segments$putative_chimeric == FALSE]))
    n_chimeric = length(unique(blast.res.segments$read_id[blast.res.segments$putative_chimeric == TRUE]))

    message(paste0(n_reads, " reads with a detected pcr adaptor"))
    message(paste0(n_non_chimeric, " non chimeric reads (", round(100*n_non_chimeric/n_reads, 2), "%)"))
    message(
      paste0(n_chimeric,
             " chimeric reads (", round(100*n_chimeric/n_reads, 2), "%) have been split into ",
             length(blast.res.segments$read_id[blast.res.segments$putative_chimeric == TRUE]), " segments"))
  }

  # 7. Rewrite the FASTQ files with the selected split-segments
  if(!silent){message("5. Create new FASTQ files containing the regular reads and splitted chimeric segments")}

  #split_and_write_fastq(FASTQ = FASTQ,
  #                      CHIMERIC_SEGMENTS_PATH = paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_chimericSegments.info.txt"),
  #                      CHIMERIC_READ_IDs = paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_chimericReads.ids.txt"),
  #                      OUTPUT_PATH = OUTPUT_DIR,
  #                      OUTPUT_SUFFIX = OUTPUT_SUFFIX)

  system(intern = FALSE, ignore.stderr = TRUE, command = paste0("seqkit grep -f ",
                                                                paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_chimericReads.ids.txt"), " ", FASTQ, " -o ", OUTPUT_DIR,
                                                                "/", OUTPUT_SUFFIX, ".chimeric_reads.fq.gz"))

  system(intern = FALSE, ignore.stderr = TRUE, command = paste0("seqkit grep -vf ",
                                                                paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_chimericReads.ids.txt"), " ", FASTQ, " -o ", OUTPUT_DIR,
                                                                "/", OUTPUT_SUFFIX, ".non_chimeric_reads.fq.gz"))

  python_script <- system.file("extdata", "split_and_write_fastq.py", package="FSnanoporeR")
  system(intern = FALSE, command = paste0(
    "python ", python_script,
    " --input_fastq ", paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, ".chimeric_reads.fq.gz"),
    " --segment_file ", paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_chimericSegments.info.txt"),
    " --output_fastq_chimeric ", paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, ".chimeric_segments.fq.gz"))
  )

  if(!silent){message("6. Finished Spliting Chimeric Reads")}

  write.table(blast.res.segments,
              file = paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_splitChimeric_reads.report.txt"),
              quote = FALSE,
              sep = "\t",
              row.names = FALSE)

  return("Finished")
}


#' Compare chimeric segments and reads using alignment data
#'
#' @description
#' Summary: Validate the split of chimeric reads into sub-segments using mapping information. Exclude in-silico chimeric.
#'
#' PART-1: Ligation Chimerics
#'
#' 1. Map chimeric segments and chimeric reads using Minimap2 ('run_minimap2()').
#'
#' 2. Extract the mapped segments (M: mapped, S: soft-clipped, H: Hard-clipped) in both primary and supplementary mappings.
#'
#' 3. Convert the read/segment mapped into GRanges objects.
#'
#' 4. For each read, compare the mapped segments to the read mapping. Assumption: the chimeric segments align at the same position as the read.
#'
#' 5. Calculate the overlap between the mapped segments and the union of the mapped segment/read. Preferred over a simple % of intersection to avoid this case:
#'
#'
#' Segment:    |===|
#'
#'
#' Read:   |===========|
#'
#'
#' Intersection shows 100% but union <50%.
#'
#' 6. For each segment / read overlap, select the one with the highest overlap.
#'
#' 7. Segment which have < "OVERLAP_THRESHOLD" with the read mapping portions are discarded.
#'
#' This severs as a validation of the chimeric split. Will exclude both false negative and positive chimeric events.
#'
#' PART-2: In-silico Chimerics
#'
#' 8. Detect segments of a same read mapping to the same genomic position (overlap >"OVERLAP_THRESHOLD").
#'
#' 9. If they display opposite orientations = in-silico chimerics.
#'
#' 10. For each in-silico pairs, select the longest segment length as the correct one.
#'
#' Final:
#'
#' Report the segment_ids which do not display enough overlap with the read mapping.
#'
#' Report the shortest segments among in-silico chimeric reads.
#'
#' This function do not discard any segment / read. It will flag them as problematic but it is up-to the user to decide if they should be removed.
#' I recommend using Dorado + duplex-tools to split and use in-silico chimerics for duplex. This should minimise the in-silico chimeric left-overs.
#'
#' @param BAM_path character. "/path/to/bam.sorted.bam" BAM must have been produced with both raw reads and chimeric splitted reads.
#' @param OUTPUT_DIR character. "/path/to/output/"
#' @param OUTPUT_SUFFIX character.
#' @param CHIMERIC_INFO character. "/path/to/_chimericSegments.info.txt" from 'run_blast_and_split_reads()'.
#'
#' @param OVERLAP_THRESHOLD numeric. Percentage of overlap between the union of chimeric segment and read mappings to consider the segment as chimeric read split as correct.
#'
#' @param silent boolean. If TRUE will print various messages.
#' @param THREADS numeric.
#'
#' @import tidyverse
#' @import ShortRead
#' @import GenomicRanges
#' @import Biostrings
#' @import S4Vectors
#' @import IRanges
#'
#' @return
#' data.frame - _chimeric_banList.txt
#'
#' data.frame - _in_silico_chimerics.log.txt
#'
#' data.frame - _mapping_chimerics.log.txt
#'
#' @export
chimeric_mapping_validation <- function(BAM_path = NULL,
                                        CHIMERIC_INFO = NULL,
                                        OUTPUT_DIR = NULL,
                                        OUTPUT_SUFFIX = NULL,
                                        THREADS = 1,
                                        OVERLAP_THRESHOLD = 0.8,
                                        silent = FALSE){


  silent = as.logical(toupper(silent))

  # 2. Load the chimeric blast infos
  if(!silent){message("1. Load Chimeric infos")}
  if(!is.data.frame(CHIMERIC_INFO)){
    chimeric.info <- read_tsv(CHIMERIC_INFO, show_col_types = FALSE) %>%
      mutate(segment_range = paste0(segment_start, "-", segment_end))
  } else {
    chimeric.info <- CHIMERIC_INFO %>%
      mutate(segment_range = paste0(segment_start, "-", segment_end))
  }

  # 3. Load the mapping information
  if(!silent){message("2. Load Mapping Informations and extract mapping ranges")}
  bam <- load_bam(BAM_path = BAM_path,
                  tidy_chimeric = TRUE,
                  clip_flag_width = 300) %>%
    filter(mapq > 20)

  if(nrow(bam) == 0){stop("No mapped read detected!")}

  # 4. Transform the mapped reads (bam) and mapped segments (blast) into GRanges objects
  # Only select the mapped portions (S/H-clipping excluded)
  if(!silent){message("3. Export to GRanges objects")}
  bam_reads.gr <- bam %>%
    # Select chimeric reads and only mapped parts
    dplyr::filter(!grepl("_chim", read_id), clip_type == "M") %>%
    mutate(segment_range = paste0(read_start_corrected, "-", read_end_corrected),
           chrom = seqnames,
           seqnames = paste0(simplified_read_id, "_", seqnames)) %>%
    select(seqnames, map_start, map_end, read_orientation, segment_range, chrom) %>%
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "seqnames", start.field = "map_start", end.field = "map_end", strand.field = "read_orientation", keep.extra.columns = TRUE)

  bam.segms.gr <- bam %>%
    # Select chimeric segments and only mapped parts
    dplyr::filter(grepl("_chim", read_id),
                  clip_type == "M",
                  sam_flag == "PRIMARY") %>%
    mutate(segment_id = simplified_read_id,
           chrom = seqnames,
           seqnames = paste0(str_replace(simplified_read_id, "_chim.*$", ""), "_", seqnames),
           segment_gr = paste0(str_replace(simplified_read_id, "_chim.*$", ""), "_", chrom, ":", map_start, "-", map_end, ":", read_orientation)) %>%
    left_join(select(chimeric.info, fastq_id_new, segment_range), by = c( "read_id" = "fastq_id_new")) %>%
    select(seqnames, map_start, map_end, read_orientation, segment_id, segment_range, chrom, segment_gr, read_id, read_group, read_length) %>%
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "seqnames", start.field = "map_start", end.field = "map_end", strand.field = "read_orientation", keep.extra.columns = TRUE)

  # 5. Compare the mapping of the chimeric reads and chimeric segments.
  # In proper read splits an overlap read/segments is expected.
  if(!silent){message("4. Compare chimeric reads/segments mappings")}

  hits.index <- suppressWarnings(GenomicRanges::findOverlaps(query = bam.segms.gr, subject = bam_reads.gr, type = "any"))
  segms.gr.hits <- bam.segms.gr[S4Vectors::queryHits(hits.index)]
  reads.gr.hits <- bam_reads.gr[S4Vectors::subjectHits(hits.index)]

  # Calculate the intersection of between the features (relative to the genome).
  # pintersect: proportion of query in subject
  p.intersect <- suppressWarnings(GenomicRanges::pintersect(segms.gr.hits, reads.gr.hits))
  pOverlap_mapping <- suppressWarnings(ShortRead::width(p.intersect) / ShortRead::width(segms.gr.hits))
  # union: portion of query in the union of the query/subject
  # Used to avoid cases where the query is smaller than the subject and fully encompassed in it.
  union_width <- suppressWarnings(GenomicRanges::punion(segms.gr.hits, reads.gr.hits))
  coverage <- ShortRead::width(p.intersect) / ShortRead::width(union_width)

  # 6. Tidy
  if(!silent){message("5. Combine and Tidy")}
  read_segm.overlaps <- data.frame(
    segment_id =  segms.gr.hits$segment_id,
    # Position in the genome
    seg.genome = paste0(segms.gr.hits$chrom, ":", GenomicRanges::ranges(segms.gr.hits), ":", ShortRead::strand(segms.gr.hits)),
    read.genome =  paste0(reads.gr.hits$chrom, ":", GenomicRanges::ranges(reads.gr.hits), ":", ShortRead::strand(reads.gr.hits)),
    # Overlaps (relative to the genome mapping)
    pintersect.genome = pOverlap_mapping,
    cov.genome = coverage,
    # Position relative to the read/segment usage
    seg.range.inRead = segms.gr.hits$segment_range,
    seg.width = GenomicRanges::width(IRanges::IRanges(segms.gr.hits$segment_range)),
    seg.length = segms.gr.hits$read_length,
    read.range = reads.gr.hits$segment_range,
    read.width = GenomicRanges::width(IRanges::IRanges(reads.gr.hits$segment_range)),
    read_id = segms.gr.hits$read_id
  )

  # 7. Select the best match between reads/segments
  # In some rare cases one chimeric segment / read can map to several places.
  # Often this is due to a misalignment of minimap2 resulting from a in-silico read not yet split.
  # As long as they have one good match segment / read I keep them.
  if(!silent){message("6. Flag disconnected chimeric read/segment mappings.")}
  read_segm.overlaps.best <- mutate(read_segm.overlaps, SUSPICIOUS_MAPPING = cov.genome < OVERLAP_THRESHOLD)

  #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

  if(!silent){message("7. Detect in-silico chimeric reads")}
  # To exclude the self-folded reads / in-silico 1d^2 chimeric
  # https://f1000research.com/articles/6-631/v2
  # https://github.com/nanoporetech/dorado
  # Dorado + duplex split_read function can be used on these cases but guppy does ignore in-silico chimeric.
  # Typically quite low but can be problematic if overloading a flowcell.

  bam.segms.gr <- bam %>%
    # Filter out problematic segments from the read/segment validation
    dplyr::filter(!read_id %in% read_segm.overlaps.best$read_id[read_segm.overlaps.best$SUSPICIOUS_MAPPING]) %>%
    # Select chimeric segments and only mapped parts
    dplyr::filter(grepl("_chim", read_id),
                  clip_type == "M") %>%
    mutate(segment_id = simplified_read_id,
           chrom = seqnames,
           seqnames = paste0(str_replace(simplified_read_id, "_chim.*$", ""), "_", seqnames),
           segment_gr = paste0(str_replace(simplified_read_id, "_chim.*$", ""), "_", chrom, ":", map_start, "-", map_end, ":", read_orientation)) %>%
    left_join(select(chimeric.info, fastq_id_new, segment_range), by = c( "read_id" = "fastq_id_new")) %>%
    select(seqnames, map_start, map_end, read_orientation, segment_id, segment_range, chrom, segment_gr, read_id, read_group, read_length) %>%
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "seqnames", start.field = "map_start", end.field = "map_end", strand.field = "read_orientation", keep.extra.columns = TRUE)

  # 8. Get the overlap between read mapping INSIDE a chimeric segment.
  # Append the read_id_chr:start-end:strand as GRanges. Use the common read_id (not unique to segment) to loop through all the segments of a read.
  # Compare mapping positions only between the same reads.
  segment.self.index <- GenomicRanges::findOverlaps(GenomicRanges::GRanges(bam.segms.gr$segment_gr), type = "any", ignore.strand = TRUE)

  self.query.index <- S4Vectors::queryHits(segment.self.index)
  self.subject.index <- S4Vectors::subjectHits(segment.self.index)
  self.gr_query.hits <- bam.segms.gr[self.query.index]
  self.gr_subject.hits <- bam.segms.gr[self.subject.index]

  # Case #1: Internal to segment ==> Great/Perfect genome match but opposite strands.
  # These cases are not dangerous, they will just map as a supplementary alignment.
  # Case #2: Between-segments ==> Great/Perfect genome match but opposite strands.
  # These cases are dangerous. This means that a read can be double-counted as it has been split into two segments.
  segment_self.overlaps <- tibble(
    read_id = self.gr_query.hits$read_id,
    simplified_read_id = str_extract(as.character(GenomicRanges::seqnames(self.gr_query.hits)), "mapread_[^_]+"),
    # Which segment # is associated to the query or subject
    query_seg_id = bam.segms.gr$read_group[self.query.index],
    subject_seg_id = bam.segms.gr$read_group[self.subject.index],
    # Chromosome
    seqnames = str_extract(as.character(GenomicRanges::seqnames(self.gr_query.hits)), "[^_]+$"),
    # Genomic position
    genome.query.range = paste0(GenomicRanges::ranges(self.gr_query.hits)),
    genome.query.strand = as.character(BiocGenerics::strand(self.gr_query.hits)),
    genome.subject.range = paste0(GenomicRanges::ranges(self.gr_subject.hits)),
    genome.subject.strand = as.character(BiocGenerics::strand(self.gr_subject.hits)),
    # Percentage of overlap between the two genomic positions
    cov.segments = BiocGenerics::width(GenomicRanges::pintersect(self.gr_query.hits, self.gr_subject.hits, ignore.strand = TRUE)) / BiocGenerics::width(GenomicRanges::punion(self.gr_query.hits, self.gr_subject.hits, ignore.strand = TRUE)),
    # Position in the chimeric read and chimeric segment width
    segm.range.inRead.query = bam.segms.gr$segment_range[self.query.index],
    segm.range.inRead.subject = bam.segms.gr$segment_range[self.subject.index],
    segm.width.query = BiocGenerics::width(IRanges::IRanges(bam.segms.gr$segment_range[self.query.index])),
    segm.width.subject = BiocGenerics::width(IRanges::IRanges(bam.segms.gr$segment_range[self.subject.index]))
  )

  in_silico_chimer <- dplyr::filter(
    segment_self.overlaps, cov.segments > 0 & query_seg_id != subject_seg_id & genome.query.strand != genome.subject.strand) %>%
    group_by(simplified_read_id) %>%
    mutate(shortest_read_id = paste0(str_extract(unique(read_id), "[^_]+"), "_chimeric=TRUE_", ifelse(segm.width.query > segm.width.subject, query_seg_id, subject_seg_id))) %>%
    ungroup()

  if(!silent){message("7. Create a ban list of suspicious chimeric reads without proper alignments")}

  banList <- tibble(read_id = character(), TYPE = character())

  # Chimeric segments without any match segment-reads
  no_match <- bam.segms.gr$read_id[-unique(c(S4Vectors::queryHits(hits.index), S4Vectors::subjectHits(hits.index)))]
  if(length(no_match) > 0){
    banList <- bind_rows(banList, data.frame(read_id = no_match, TYPE = "NO_MATCH"))
  }

  # Chimeric segments with match segment-reads below threshold
  partial_match <- unique(read_segm.overlaps.best$read_id[read_segm.overlaps.best$cov.genome < OVERLAP_THRESHOLD])
  if(length(partial_match) > 0){
    banList <- bind_rows(banList, data.frame(read_id = partial_match, TYPE = "PARTIAL_MATCH"))
  }

  # In-silico chimeric
  in_silico <- unique(in_silico_chimer$shortest_read_id)
  if(length(in_silico) > 0){
    banList <- bind_rows(banList, data.frame(read_id = in_silico, TYPE = "IN_SILICO")) %>%
      mutate(TYPE = ifelse(read_id %in% read_id[TYPE == "IN_SILICO"], "IN_SILICO", TYPE))

  }

  message(paste0(length(unique(banList$read_id)),
                 " chimeric segments detected with abnormal mappings - exclude them! (",
                 round(100*length(unique(banList$read_id))/length(unique(bam.segms.gr$read_id)), 2), "% of the chimeric segments)"))
  message(paste0("Among which ", length(unique(banList$read_id[banList$TYPE == "IN_SILICO"])),
                 " in-silico chimeric read split (",
                 round(100*length(unique(banList$read_id[banList$TYPE == "IN_SILICO"]))/length(unique(bam.segms.gr$segment_id)), 2), "% of the chimeric segments)"))

  if(!silent){message("8. Save Results")}
  write.table(banList,
              file = paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_chimeric_banList.txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

  write.table(read_segm.overlaps,
              file = paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_mapping_chimerics.log.txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

  write.table(segment_self.overlaps,
              file = paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_in_silico_chimerics.log.txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

  return("Finished")

}
