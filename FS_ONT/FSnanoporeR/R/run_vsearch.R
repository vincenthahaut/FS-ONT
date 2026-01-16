#' run_vsearch but on the read ends
#'
#' @description
#' Search for a specific set of barcode sequences in FASTQ reads. Vsearch is preferred over BLAST here because of the way it handles Ns and the output format.
#' Altogether this provides a slightly better performance than BLAST at retrieving the primers.
#' The difference between run_vsearch and run_vsearch2 is that the later first extract the WINDOW bp (ex 200bp) from the read end/start and then look for the barcode in that specific region while the other does it accross the entire read body.
#' Probably feasible with BLAST but quite simple with VSEARCH once chimeric reads have been split.
#' Warning: Will only find one hit per read which prevents it from detecting chimeric molecules other than those with two different barcodes!
#'
#' Requires vsearch installed and in your path (https://github.com/torognes/vsearch).
#'
#' @param FASTQ character. /path/to/the/input/fastq.gz. Gziped (.gz) or not (detected with grepl).
#' @param OUTPUT_DIR character. /path/to/the/output_directory/
#' @param FORWARD_SEQ character. Sequence of the forward sequence to search for (5==>3).
#' @param REVERSE_SEQ character. Sequence of the reverse sequence to search for (5==>3).
#' @param OUTPUT_SUFFIX character. Sample id to label the outputs.
#' @param SEQ_IDENTITY numeric. Minimal identity required between the barcode sequence and the query.
#' @param MIN_SEQ_LENGTH numeric. Minimal match length between the barcode sequence and the query.
#' @param WINDOW numeric. Run vsearch only on the read start/end windows.
#' @param silent boolean. Display messages or not.
#'
#' @param THREADS numeric. Number of threads to run vsearch
#'
#' @details Vsearch fields: https://www.drive5.com/usearch/manual/userfields.html
#'
#' @returns Write the results to paste0(OUTPUT_DIR, "vsearch_", SAMPLE_ID, ".txt").
#'
#' @export
run_vsearch <- function(FASTQ = NULL,
                        OUTPUT_DIR = "./",
                        FORWARD_SEQ = NULL,
                        REVERSE_SEQ = NULL,
                        SEQ_IDENTITY = 0.75,
                        OUTPUT_SUFFIX = "bc" ,
                        MIN_SEQ_LENGTH = 30,
                        WINDOW = 200,
                        THREADS = 1,
                        silent = FALSE){
  
  if(grepl("not found", system(intern = TRUE, "type vsearch &>/dev/null"))){
    stop("Please install vsearch")
  }
  if(grepl("not found", system(intern = TRUE, "type seqkit &>/dev/null"))){
    stop("Please install seqkit")
  }
  
  # 1. Pre-filter the by read length and convert to fasta
  # Not per-se mandatory but prevents some unwanted behaviour when barcode detection / read trimming
  # Throws a warning "[WARN] you may switch on flag -g/--remove-gaps to remove spaces" - Normal and should be left.
  invisible(system(intern = FALSE, ignore.stdout = silent, ignore.stderr = silent, paste0(
    "zcat ", paste0(FASTQ, collapse = " "), " | ",
    "seqkit seq --threads ", THREADS, " --min-len ", WINDOW, " - | ",
    "seqkit fq2fa --threads ", THREADS, " - -o ", OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_vsearch_filtered.fa")))
  
  # 2. Extract the read start / ends
  system(intern = FALSE, ignore.stderr = silent,
         command = paste0("seqkit subseq --threads ", THREADS, " -r 1:", WINDOW, " ", OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_vsearch_filtered.fa > ",
                          OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_vsearch_readstart.fa")
  )
  
  system(intern = FALSE, ignore.stderr = silent,
         command = paste0("seqkit subseq --threads ", THREADS, " -r -", WINDOW ,":-1 ", OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_vsearch_filtered.fa > ",
                          OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_vsearch_readend.fa")
  )
  
  # 3. Create the primer database
  vsearch_db_fwd <- paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "PRIMER_vsearch_fwd_db.fa")
  write(">BC1_fwd", vsearch_db_fwd)
  write(REVERSE_SEQ, vsearch_db_fwd, append = TRUE)
  
  vsearch_db_rev <- paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "PRIMER_vsearch_rev_db.fa")
  write(">BC1_rev", vsearch_db_rev)
  write(as.character(ShortRead::reverseComplement(Biostrings::DNAString(FORWARD_SEQ))), vsearch_db_rev, append = TRUE)
  
  # 4. RUN VSEARCH
  system(intern = FALSE, ignore.stdout = silent, ignore.stderr = silent,
         command = paste(
           "vsearch --usearch_global", paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_vsearch_readstart.fa"),
           "--threads", THREADS,
           #"--minseqlength", MIN_SEQ_LENGTH,
           "--maxaccepts 5",
           "--strand plus",
           "--wordlength 3",
           "--minwordmatches 12",
           "--userfields 'query+target+id+alnlen+mism+opens+qilo+qihi+qstrand+tilo+tihi+ql+tl+qrow+trow'",
           "--mincols", MIN_SEQ_LENGTH,
           "--userout", paste0(OUTPUT_DIR, "/vsearch_readstart_", OUTPUT_SUFFIX, ".txt"),
           "--db", vsearch_db_fwd,
           "--id", SEQ_IDENTITY,
           ifelse(silent, "--quiet", "")
         )
  )
  
  system(intern = FALSE, ignore.stdout = silent, ignore.stderr = silent,
         command = paste(
           "vsearch --usearch_global", paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_vsearch_readend.fa"),
           "--threads", THREADS,
           #"--minseqlength", MIN_SEQ_LENGTH,
           "--maxaccepts 5",
           "--strand plus",
           "--wordlength 3",
           "--minwordmatches 12",
           "--userfields 'query+target+id+alnlen+mism+opens+qilo+qihi+qstrand+tilo+tihi+ql+tl+qrow+trow'",
           "--mincols 9",
           "--userout", paste0(OUTPUT_DIR, "/vsearch_readend_", OUTPUT_SUFFIX, ".txt"),
           "--db", vsearch_db_rev,
           "--id", SEQ_IDENTITY,
           ifelse(silent, "--quiet", "")
         )
  )
  
  # 5. Combine results
  # Would be nice to do one liner but R did not like <()
  file1 <- paste0(OUTPUT_DIR, "/vsearch_readstart_", OUTPUT_SUFFIX, ".txt")
  file2 <- paste0(OUTPUT_DIR, "/vsearch_readend_", OUTPUT_SUFFIX, ".txt")
  output <- paste0(OUTPUT_DIR, "/vsearch_", OUTPUT_SUFFIX, ".txt")
  tmp_file1 <- tempfile(fileext = ".txt")
  tmp_file2 <- tempfile(fileext = ".txt")
  
  system(sprintf("awk '{print $0 \"\\tSTART\"}' %s > %s", file1, tmp_file1))
  system(sprintf("awk 'NR!=1{print $0 \"\\tEND\"}' %s > %s", file2, tmp_file2))
  
  system(sprintf("cat %s %s > %s", tmp_file1, tmp_file2, output))
  
  file.remove(tmp_file1)
  file.remove(tmp_file2)
  
  # 6. Add column names
  system(intern = FALSE,
         command = paste(
           "sed -i",
           "'1i read_id	barcode	identity	bc_match_l	mism	gap_open	read_start	read_end	read_strand	bc_start	bc_end	query_l	bc_l	query_seq_aln	bc_seq_aln	read_section'",
           paste0(OUTPUT_DIR, "/vsearch_", OUTPUT_SUFFIX, ".txt")
         )
  )
  
  # 7. Clean-up
  system(intern = FALSE,
         command = paste("rm",
                         paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_vsearch_filtered.fa"),
                         paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_vsearch_filtered.fa.seqkit.fai"),
                         paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "PRIMER_vsearch_fwd_db.fa"),
                         paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "PRIMER_vsearch_rev_db.fa"),
                         paste0(OUTPUT_DIR, "/vsearch_readstart_", OUTPUT_SUFFIX, ".txt"),
                         paste0(OUTPUT_DIR, "/vsearch_readend_", OUTPUT_SUFFIX, ".txt"), "&")
  )
  
}