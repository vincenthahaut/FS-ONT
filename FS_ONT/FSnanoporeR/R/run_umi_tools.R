#' run_umi_tools
#'
#' Run umi_tools on one or multiple bam files using first the assigned gene tag and then the transcript gene tag from isoquant.
#'
#' @param BAM character. "/path/to/the/file_isoquant_umi.bam" converted using convert_isoquant_UMItools.
#' @param OUTPUT_SUFFIX character. Sample id to label the outputs.
#' @param OUTPUT_DIR character. "/path/to/the/output_directory"
#' @param UMI_TOOLS_BIN, character. "Path/to/umi_tools" binaries. If not set, default to "umi_tools".
#'
#' @export
run_umi_tools <- function(BAM = NULL, OUTPUT_SUFFIX = NULL, OUTPUT_DIR = NULL, silent = TRUE, UMI_TOOLS_BIN = "umi_tools"){

  if(!file.exists(OUTPUT_DIR)){dir.create(OUTPUT_DIR)}

  # Run UMI Tools using the gene tag (XG)
  system(intern=FALSE, command = paste0(
    UMI_TOOLS_BIN, " count --per-gene --gene-tag=XG --assigned-status-tag=GS --per-cell --cell-tag=CB --extract-umi-method=tag --umi-tag=UB ",
    "--log=", paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_isoquant_umi_gene_logs.tsv"),
    " -I ", BAM, " -S ",
    paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_isoquant_umi_gene_counts.tsv")
    )
  )

  # Run UMI Tools using the transcript tag (XT)
  system(intern=FALSE, command = paste0(
    UMI_TOOLS_BIN, " count --per-gene --gene-tag=XT --assigned-status-tag=TS --per-cell --cell-tag=CB --extract-umi-method=tag --umi-tag=UB ",
    "--log=", paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_isoquant_umi_transcript_logs.stringent.unambiguous.tsv"),
    " -I ", BAM, " -S ",
    paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_isoquant_umi_transcript_counts.stringent.unambiguous.tsv")
    )
  )

  # Run UMI Tools using the transcript tag (XI)
  system(intern=FALSE, command = paste0(
    UMI_TOOLS_BIN, " count --per-gene --gene-tag=XI --assigned-status-tag=IS --per-cell --cell-tag=CB --extract-umi-method=tag --umi-tag=UB ",
    "--log=", paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_isoquant_umi_transcript_logs.stringent.tsv"),
    " -I ", BAM, " -S ",
    paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_isoquant_umi_transcript_counts.stringent.tsv")
    )
  )

  # Run UMI Tools using gene tag (RS)
  system(intern=FALSE, command = paste0(
    UMI_TOOLS_BIN, " count --per-gene --gene-tag=XR --assigned-status-tag=RS --per-cell --cell-tag=CB --extract-umi-method=tag --umi-tag=UB ",
    "--log=", paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_isoquant_umi_gene_logs.stringent.tsv"),
    " -I ", BAM, " -S ",
    paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_isoquant_umi_gene_counts.stringent.tsv")
  )
  )

}

#' tidy_umi_tools
#'
#' Create a matrix of counts using umi_Tools results
#'
#' @param OUTPUT_SUFFIX character. Sample id to label the outputs.
#' @param OUTPUT_DIR character. "/path/to/the/output_directory"
#' @param UMI_TOOLS_RESULTS_FOLDER character. Sample id to label the outputs.
#'
#' @return OUTPUT_DIR/OUTPUT_SUFFIX_isoquant_umitools_counts.gene.tsv.gz and OUTPUT_DIR/OUTPUT_SUFFIX_isoquant_umitools_counts.transcripts.tsv.gz files
#'
#' @export
tidy_umi_tools <- function(UMI_TOOLS_RESULTS_FOLDER = NULL, OUTPUT_DIR = NULL, OUTPUT_SUFFIX = NULL){

  # Genes
  gene.df <- lapply(list.files(UMI_TOOLS_RESULTS_FOLDER, full.names = TRUE, pattern = "_isoquant_umi_gene_counts.tsv"), function(x) read_tsv(x, show_col_types = FALSE))
  gene.df <- gene.df[sapply(gene.df, nrow)  >0 ]

  if(length(gene.df) > 0){
    cell_ids <- as.vector(unlist(lapply(gene.df, function(x) unique(x["cell"]))))

    gene.mat <- suppressMessages(lapply(gene.df, function(x) select(x, - cell)) %>% purrr::reduce(full_join, by = "gene"))
    colnames(gene.mat) <- c("gene", cell_ids)

    gene_ids <- gene.mat$gene
    gene_counts <- select(gene.mat, -gene) %>%
      mutate_all(~replace_na(.,0))

    gene.mat <- bind_cols(data.frame(gene = gene_ids), gene_counts)

    out.gene <- paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_isoquant_umitools_counts.gene.tsv")
    if(file.exists(out.gene)){file.remove(out.gene)}

    write_tsv(gene.mat, out.gene, col_names = TRUE)

  }

  # Genes - stringent
  gene.df <- lapply(list.files(UMI_TOOLS_RESULTS_FOLDER, full.names = TRUE, pattern = "_isoquant_umi_gene_counts.stringent.tsv"), function(x) read_tsv(x, show_col_types = FALSE))
  gene.df <- gene.df[sapply(gene.df, nrow)  >0 ]

  if(length(gene.df) > 0){
    cell_ids <- as.vector(unlist(lapply(gene.df, function(x) unique(x["cell"]))))

    gene.mat <- suppressMessages(lapply(gene.df, function(x) select(x, - cell)) %>% purrr::reduce(full_join, by = "gene"))
    colnames(gene.mat) <- c("gene", cell_ids)

    gene_ids <- gene.mat$gene
    gene_counts <- select(gene.mat, -gene) %>%
      mutate_all(~replace_na(.,0))

    gene.mat <- bind_cols(data.frame(gene = gene_ids), gene_counts)

    out.gene <- paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_isoquant_umitools_counts.gene.stringent.tsv")
    if(file.exists(out.gene)){file.remove(out.gene)}

    write_tsv(gene.mat, out.gene, col_names = TRUE)

  }

  # Transcripts (unique only)
  transcript.df <- lapply(list.files(UMI_TOOLS_RESULTS_FOLDER, full.names = TRUE, pattern = "_isoquant_umi_transcript_counts.stringent.unambiguous.tsv"), function(x) read_tsv(x, show_col_types = FALSE))

  transcript.df <- transcript.df[sapply(transcript.df, nrow)>0]

  if(length(transcript.df) > 0){

    cell_ids <- as.vector(unlist(lapply(transcript.df, function(x) unique(x["cell"]))))
    transcript.mat <- suppressMessages(lapply(transcript.df, function(x) select(x, - cell)) %>% purrr::reduce(full_join, by = "gene"))
    colnames(transcript.mat) <- c("gene", cell_ids)

    transcript.mat_ids <- transcript.mat$gene
    transcript.mat_counts <- select(transcript.mat, -gene) %>%
      mutate_all(~replace_na(.,0))

    transcript.mat <- bind_cols(data.frame(gene = transcript.mat_ids), transcript.mat_counts)

    # Save
    out.transcript <- paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_isoquant_umitools_counts.stringent.unambiguous.tsv")

    if(file.exists(out.transcript)){file.remove(out.transcript)}

    write_tsv(transcript.mat, out.transcript, col_names = TRUE)
  }

  # Transcripts (unique and inconsistent)
  transcript.df <- lapply(list.files(UMI_TOOLS_RESULTS_FOLDER, full.names = TRUE, pattern = "_isoquant_umi_transcript_counts.stringent.tsv"), function(x) read_tsv(x, show_col_types = FALSE))

  transcript.df <- transcript.df[sapply(transcript.df, nrow)>0]

  if(length(transcript.df) > 0){

    cell_ids <- as.vector(unlist(lapply(transcript.df, function(x) unique(x["cell"]))))
    transcript.mat <- suppressMessages(lapply(transcript.df, function(x) select(x, - cell)) %>% purrr::reduce(full_join, by = "gene"))
    colnames(transcript.mat) <- c("gene", cell_ids)

    transcript.mat_ids <- transcript.mat$gene
    transcript.mat_counts <- select(transcript.mat, -gene) %>%
      mutate_all(~replace_na(.,0))

    transcript.mat <- bind_cols(data.frame(gene = transcript.mat_ids), transcript.mat_counts)

    # Save
    out.transcript <- paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_isoquant_umitools_counts.transcripts.stringent.tsv")

    if(file.exists(out.transcript)){file.remove(out.transcript)}

    write_tsv(transcript.mat, out.transcript, col_names = TRUE)
  }
}


#' update_read_assignments
#'
#' Update the read assignment file from Isoquant to add new SAM tags
#'
#' NB: There is definitely a faster way to do it - I will explore options later.
#'
#' @param ISOQUANT_READ_ASSIGNMENTS character. "/path/to/the/input/*.read_assignments.tsv" from isoquant.
#' @param OUTPUT_SUFFIX character. Sample id to label the outputs.
#' @param OUTPUT_DIR character. "/path/to/the/output/"
#' @param silent boolean. Display messages or not.
#' @param mode.bam "NONE", "MONOMER", or "TRIMER"
#'
#' @export
update_read_assignments <- function(ISOQUANT_READ_ASSIGNMENTS = NULL, mode.bam = NULL, OUTPUT_SUFFIX = NULL, OUTPUT_DIR = NULL, silent = FALSE){

  if(!silent){message("1. Update Read Tag File")}

  iso_reads <- vroom::vroom(ISOQUANT_READ_ASSIGNMENTS, comment = "#", col_names = FALSE, show_col_types = FALSE)
  colnames(iso_reads) <- c("read_id", "chr", "strand", "isoform_id", "gene_id", "assignment_type", "assignment_events", "exons", "additional_info")

  iso_reads$classification <- str_extract(iso_reads$additional_info, "Classification=.*[^;]") %>%
    str_replace("Classification=", "")
  iso_reads$canonical <- str_extract(iso_reads$additional_info, "(?<=Canonical=)[^;]+")

  if(!silent){message("2. Update Inconsistent Reads")}
  # Significant Inconsistencies
  major_inconsistencies <- c(
    # Introns
    "intron_migration",
    "intron_alternation",
    # Exons
    "major_exon_elongation",
    "mutually_exclusive_exons",
    "exon_skipping",
    "exon_merge",
    "exon_gain",
    "exon_detach",
    # Splicing
    "alt_donor_site_novel",
    "alt_acceptor_site_novel",
    # Other
    "alternative_structure",
    # Typical from highly truncated molecules / gDNA
    "mono_exonic",
    "antisense",
    "incomplete_intron_retention",
    "internal_polya"
  )

  # Filter unique / inconsistent reads that do not match these major inconsistencies

  # A. Unique reads with major mapping inconsistencies
  iso_reads <- mutate(iso_reads, index = 1:n())
  iso_reads.filtered.unique <- dplyr::filter(iso_reads,
                  assignment_type %in% c("unique","unique_minor_difference"),
                  classification != "mono_exonic_match",
                  grepl(paste0(major_inconsistencies, collapse = "|"), assignment_events)
  )

  iso_reads$assignment_type[iso_reads$index %in% iso_reads.filtered.unique$index] <- "unique_banned"

  # B. Inconsistent with no major mapping alterations
  iso_reads.filtered.inconsistent <- dplyr::filter(iso_reads,
                                      assignment_type == "inconsistent",
                                      !grepl(paste0(major_inconsistencies, collapse = "|"), assignment_events)
  )

  iso_reads$assignment_type[iso_reads$index %in% iso_reads.filtered.inconsistent$index] <- "inconsistent_recovered"

  # C. Any reads fully unspliced but not belonging to a mono_exonic_match
  iso_reads <- mutate(iso_reads, index = 1:n())
  iso_reads.filtered.unspliced <- dplyr::filter(iso_reads,
                                                assignment_type %in% c("unique","unique_minor_difference", "inconsistent"),
                                                classification != "mono_exonic_match",
                                                canonical == "Unspliced"
  )

  iso_reads$assignment_type[iso_reads$index %in% iso_reads.filtered.unspliced$index] <- "unspliced"

  iso_reads.filtered.inconsistent <- iso_reads.filtered.unspliced <- iso_reads.filtered.unique <- NULL
  # NB: no collision between unique and inconsistent possible.

  if(!silent){message("3. Create new SAM Tags - aggregate duplicates")}
  isoform_id <- dplyr::filter(iso_reads, assignment_type %in% c("unique", "unique_minor_difference")) %>%
    #dplyr::filter(canonical == "True")
    dplyr::select(read_id, isoform_id) %>%
    distinct() %>%
    # Deal with reads assigned to multiple transcripts
    dplyr::filter(!read_id %in% read_id[duplicated(read_id)]) %>%
    dplyr::rename("XT" = isoform_id) %>%
    mutate(TS = "Assigned")

  isoform_id_nonfull <- dplyr::filter(iso_reads, assignment_type %in% c("unique", "unique_minor_difference", "inconsistent_recovered")) %>%
    #dplyr::filter(canonical == "True")
    dplyr::select(read_id, isoform_id) %>%
    distinct() %>%
    # Deal with reads assigned to multiple transcripts
    dplyr::filter(!read_id %in% read_id[duplicated(read_id)]) %>%
    dplyr::rename("XI" = isoform_id) %>%
    mutate(IS = "Assigned")

  # Use also "inconsistent" categories for the gene ids
  gene_ids <- dplyr::filter(iso_reads, assignment_type %in% c("inconsistent", "unique", "unique_minor_difference","unique_banned", "inconsistent_recovered")) %>%
    #dplyr::filter(canonical == "True")
    dplyr::select(read_id, gene_id) %>%
    distinct() %>%
    # A few duplicated read ids may occur
    # Most of which are related to near-perfect gene overlaps in the annotation (ex: PINK1/PINK1-AS1, ASDURF/ASNSD1) or problematic reads left
    # Mask them to avoid confusions
    dplyr::filter(!read_id %in% read_id[duplicated(read_id)]) %>%
    dplyr::rename("XG" = gene_id) %>%
    mutate(GS = "Assigned")

  # Use also "inconsistent" categories for the gene ids
  gene_ids_stringent <- dplyr::filter(iso_reads, assignment_type %in% c("unique", "unique_minor_difference", "inconsistent_recovered")) %>%
    #dplyr::filter(canonical == "True")
    dplyr::select(read_id, gene_id) %>%
    distinct() %>%
    # A few duplicated read ids may occur
    # Most of which are related to near-perfect gene overlaps in the annotation (ex: PINK1/PINK1-AS1, ASDURF/ASNSD1) or problematic reads left
    # Mask them to avoid confusions
    dplyr::filter(!read_id %in% read_id[duplicated(read_id)]) %>%
    dplyr::rename("XR" = gene_id) %>%
    mutate(RS = "Assigned")

  # Create Tag file
  # XT - associated isoform id - uniquely assigned only (only relevant for UMI)
  # TS - associated isoform id - Assigned or Unassigned
  # XI - associated isoform id - full / incomplete
  # IS - assigned isoform info - Assigned or Unassigned
  # XG - associated gene id
  # GS - assigned gene info - Assigned or Unassigned
  # XR - associated gene id - stringent
  # RS - assigned gene info - Assigned or Unassigned
  # CB - cell barcode
  # UB - UMI
  # IF - Isoquant Flag(s)
  # IC - Isoquant classification
  # IN - Check canonical splicing
  if(!silent){message("4. Create read assignment file")}

  tags <- tibble(
    read_id = iso_reads$read_id,
    IF = iso_reads$assignment_type,
    IC = iso_reads$classification,
    IN = iso_reads$canonical) %>%
    distinct() %>%
    mutate(
      IF = ifelse(read_id %in% unique(read_id[duplicated(read_id)]), "mulitple", IF),
      IN = ifelse(read_id %in% unique(read_id[duplicated(read_id)]), "mulitple", IN),
      IC = ifelse(read_id %in% unique(read_id[duplicated(read_id)]), "mulitple", IC)
    ) %>%
    distinct()

  if(mode.bam %in% c("MONOMER", "TRIMER")){
    if(!silent){message("4.1. FLASH-seq-ONT UMI Mode detected- Deal with UMI")}
    tags$UB <- str_extract(pattern = "[A|T|C|G]{8}$", tags$read_id)
    tags <- dplyr::filter(tags, !is.na(UB))
  } else {
    print("no UMI detected or 10x - skip filtering by UMI")
  }

  if(!silent){message("5. Combine new read assignment file with new SAM tags - slow")}

  # Add Tags
  tags <- tags %>%
    left_join(isoform_id, by = "read_id") %>%
    left_join(isoform_id_nonfull, by = "read_id") %>%
    left_join(gene_ids, by = "read_id") %>%
    left_join(gene_ids_stringent, by = "read_id")

  # Stringent genes (not all inconsistent)
  tags$RS[is.na(tags$RS)] <- "Unassigned"
  tags$XR[is.na(tags$XR)] <- "None"
  # Stringent Transcripts (full)
  tags$TS[is.na(tags$TS)] <- "Unassigned"
  tags$XT[is.na(tags$XT)] <- "None"
  # Genes
  tags$GS[is.na(tags$GS)] <- "Unassigned"
  tags$XG[is.na(tags$XG)] <- "None"
  # transcripts (incomplete also)
  tags$IS[is.na(tags$IS)] <- "Unassigned"
  tags$XI[is.na(tags$XI)] <- "None"

  if(!silent){message("6. Save file")}
  vroom::vroom_write(tags, paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_read_assignments.updated.txt"), delim= "\t")

}

#' convert_isoquant_UMItools2
#'
#' Run a python script to add new SAM tags to a BAM. Works by crossing the BAM files with a txt file containing the new tags in columns (2 characters columns) and the read_id column (first col)
#'
#' @param ISOQUANT_BAM_FOLDER character. "/path/to/the/aux/" from isoquant (contains the BAM Files)
#' @param mode.bam "10x" or "NONE", "MONOMER", or "TRIMER"
#' @param OUTPUT_SUFFIX character. Sample id to label the outputs.
#' @param OUTPUT_DIR character. "/path/to/the/output/"
#' @param UPDATED_READ_ASSIGN character. "/path/to/updated_read_assignments.txt" from update_read_assignments()
#' @param THREADS numeric. Number of parallel sessions.
#' @param silent boolean. Suppress some of the messages.
#'
#' @export
convert_isoquant_UMItools2 <- function(ISOQUANT_BAM_FOLDER = NULL, silent = TRUE, mode.bam = "NONE", UPDATED_READ_ASSIGN = NULL, OUTPUT_SUFFIX = NULL, OUTPUT_DIR = NULL, THREADS = 4){

  if(!file.exists(OUTPUT_DIR)){dir.create(OUTPUT_DIR)}

  if(mode.bam == "10x"){
    if(!silent){message("1. Start Converting BAM - 10x Mode")}

    bams.in <- list.files(ISOQUANT_BAM_FOLDER, full.names = TRUE, pattern = "merged.bam$")
    bams.out <- paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_isoquant_umi.bam")

    # Cell tags are already present in the BAM
    python_script <- system.file("extdata", "add_tags_10x.py", package="FSnanoporeR")
    system(intern = FALSE, command = paste0(
      "python ", python_script,
      " --bam_in ", bams.in,
      " --bam_out ", bams.out,
      " --tag_file ", UPDATED_READ_ASSIGN)
    )

  } else {

    if(!silent){message("1. Start Converting BAM - FLASH-seq-ONT mode")}

    # Full-path crashes due to arg. too long
    #bams.in <- list.files(ISOQUANT_BAM_FOLDER, full.names = TRUE, pattern = ".bam$")
    setwd(ISOQUANT_BAM_FOLDER)
    bams.in <- list.files(ISOQUANT_BAM_FOLDER, pattern = ".bam$")
    bc <- str_extract(pattern = "BC\\d+", bams.in)
    # Remove chimeric / undetermined
    bams.in <- bams.in[!is.na(bc)]
    bc <- bc[!is.na(bc)]
    bams.out <- paste0(OUTPUT_DIR, OUTPUT_SUFFIX, "_", bc, "_isoquant_umi.bam")
    if(!silent){message(paste0("Found: ", length(bams.in), " files to process"))}

    python_script <- system.file("extdata", "add_tags.py", package="FSnanoporeR")
    system(intern = FALSE, command = paste0(
      "python ", python_script,
      " --threads ", THREADS,
      " --bam_in ", paste0(bams.in, collapse = " "),
      " --bam_out ", paste0(bams.out, collapse = " "),
      " --cell_barcodes ", paste0(bc, collapse = " "),
      " --tag_file ", UPDATED_READ_ASSIGN)
    )
  }
}


#' @title Isoquant Update Transcript Count Matrix
#'
#' @description
#' Using the updated read assignment file, create new grouped count matrixes for gene stringent (XR) and transcripts stringents but full/partial spliced (XI).
#'
#' @param UPDATED_READ_ASSIGN character. "/path/to/updated_read_assignments.txt" from update_read_assignments()
#' @param read_bc character. "/path/to/FSnanopore_detected_barcodes.txt" (contains columns named "read_id", "index_unique")
#' @param output_dir character. Output directory of the other functions.
#' @param output_suffix character. Sample ID.
#' @param silent boolean. TRUE or FALSE.
#'
#' @export
isoquant_count_matrix_inconsistent <- function(UPDATED_READ_ASSIGN = NULL, read_bc = "detected_barcodes.txt", OUTPUT_SUFFIX = NULL, OUTPUT_DIR = NULL, silent = TRUE){

  if(!silent){message("1. Update non-UMI count matrices - read updated read assignment file")}

  # 1. Load files and associate cell barcodes to reads
  iso_reads <- vroom::vroom(UPDATED_READ_ASSIGN, show_col_types = FALSE) %>%
    # For UMI (append to the name as _XXX)
    mutate(read_id = str_replace(read_id, "_[A|T|C|G]{8}$", ""))

  if(!silent){message("2. Add cell barcodes")}

  iso_reads <- left_join(iso_reads,
                         vroom::vroom(read_bc, show_col_types = FALSE, col_names = TRUE) %>%
                           filter(!index_unique == "Undetermined") %>%
                           select(read_id, index_unique),
                         by = c("read_id"))

  if(!silent){message("3. Create new transcript matrix")}

  tx_matrix <- iso_reads %>%
    filter(IS != "Unassigned") %>%
    select(XI, index_unique) %>%
    group_by(XI, index_unique) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    pivot_wider(values_from = n, names_from = index_unique, values_fill = 0) %>%
    dplyr::rename("transcript_id" = XI)

  write.table(tx_matrix, paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, ".transcript_grouped_counts.with_selected_inconsistent.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

  if(!silent){message("4. Create new gene matrix")}

  gene_matrix <- iso_reads %>%
    filter(RS != "Unassigned") %>%
    select(XR, index_unique) %>%
    group_by(XR, index_unique) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    pivot_wider(values_from = n, names_from = index_unique, values_fill = 0) %>%
    dplyr::rename("gene_id" = XR)

  write.table(gene_matrix, paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, ".genes_grouped_counts.with_selected_inconsistent.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

}



