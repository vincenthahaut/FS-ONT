#' @title random FASTQ sampling
#'
#' @description
#' Randomly sample a defined number of reads from a FASTQ file. Use seqkit sample
#'
#' @param fq.input character. FASTQ path (gz or not - detected with grepl).
#' @param n_reads numeric. Number of reads to downsample.
#' @param OUTPUT_DIR character. Output directory of the other functions.
#' @param OUTPUT_SUFFIX character. Sample ID.
#'
#' @export
random_sampling <- function(fq.input = NULL, n_reads = NULL, OUTPUT_SUFFIX = NULL, OUTPUT_DIR = NULL){

  system(intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE, command =
           paste0(
             ifelse(grepl(".gz$", fq.input), "zcat ", "cat "), fq.input,
             "| seqkit sample -n ", n_reads, " -o ", OUTPUT_DIR, "/", OUTPUT_SUFFIX, ".", n_reads, ".fq.gz")
  )

}


#' @title random FASTQ sampling
#'
#' @description
#' Run isoquant on the samples.
#'
#' @param downsample_n_reads numeric or "". If numeric, will resample X reads from each cell.
#' @param premapping boolean. If FALSE, let isoquant map your samples. If TRUE, map the samples in parallel prior to isoquant.
#' @param output_dir character. Output directory of the other functions.
#' @param output_suffix character. Sample ID.
#' @param isoquant_fa character. /path/to/reference/genome.fa.
#' @param isoquant_gtf character. /path/to/genome/annotation.gtf.
#' @param max_memory see run_minimap2.
#' @param minimap2_binaries see run_minimap2.
#' @param minimap2_ref see run_minimap2.
#' @param minimap2_bed see run_minimap2.
#' @param transcript_counts_type character. Accepts "unique_only" or "with_inconsistent". See isoquant.
#'
#' @export
run_isoquant <- function(downsample_n_reads = "", output_dir = NULL, output_suffix = NULL, transcript_counts_type = "unique_only", premapping = TRUE, n_cores = NULL, isoquant_fa = NULL, isoquant_gtf = NULL, max_memory = NULL, minimap2_binaries = NULL, minimap2_ref = NULL, minimap2_bed = NULL){

  # 1. Create FASTQ list for isoquant and if needed downsample reads

  # 1.1. If downsampling is required
  if(is.numeric(downsample_n_reads)){

    # 1.1.1. Get Fastq files
    fq.list <- list.files(paste0(output_dir, "/", output_suffix, "_FLASHseqONT_output/", output_suffix, "_demultiplexing_pass/"), full.names = TRUE)
    fq.ids <- list.files(paste0(output_dir, "/", output_suffix, "_FLASHseqONT_output/", output_suffix, "_demultiplexing_pass/"), full.names = FALSE)

    fq.ids <- str_replace_all(fq.ids, pattern = paste0(output_suffix, "_|_demultiplexing_pass.fq.gz"), "")

    # 1.1.2. Create Folder
    out.demultiplexing <- paste0(output_dir, "/", output_suffix, "_FLASHseqONT_output/", output_suffix, "_demultiplexing_pass_", downsample_n_reads, "reads/")
    if(file.exists(out.demultiplexing)){unlink(out.demultiplexing, recursive = TRUE)}
    dir.create(out.demultiplexing)

    # 1.1.3. Filter out samples with too few reads (seqkit doesn't do it)
    n_reads_per_bc <- read_tsv(paste0(output_dir, "/", output_suffix, "_FLASHseqONT_output/", output_suffix, "_read_length.report.txt"), show_col_types = FALSE) %>%
      group_by(cell_barcode) %>%
      summarise(n = n()) %>%
      ungroup()

    index <- fq.ids %in% n_reads_per_bc$cell_barcode[n_reads_per_bc$n > downsample_n_reads]
    fq.ids <- fq.ids[index]
    fq.list <- fq.list[index]

    # 1.1.4. Resample
    parallel::mclapply(1:length(fq.list), mc.silent = silent, mc.cores = n_cores, function(x)
      random_sampling(fq.input = fq.list[x], n_reads = downsample_n_reads, OUTPUT_SUFFIX = paste0(output_suffix, "_", fq.ids[x]), OUTPUT_DIR = out.demultiplexing)
    )

    # 1.1.5. Update fq list for isoquant
    fq.list <- list.files(out.demultiplexing, full.names = TRUE)
    fq.list.path <- paste0(output_dir, "/", output_suffix, "_FLASHseqONT_output/", output_suffix, "_fq.list")
    write.table(fq.list, fq.list.path, row.names = FALSE, quote = FALSE, col.names = FALSE)

  } else {

    # 1.2. Without downsampling
    fq.list <- list.files(paste0(output_dir, "/", output_suffix, "_FLASHseqONT_output/", output_suffix, "_demultiplexing_pass/"), full.names = TRUE)
    fq.ids <- list.files(paste0(output_dir, "/", output_suffix, "_FLASHseqONT_output/", output_suffix, "_demultiplexing_pass/"), full.names = FALSE)
    fq.list.path <- paste0(output_dir, "/", output_suffix, "_FLASHseqONT_output/", output_suffix, "_fq.list")
    fq.ids <- str_replace_all(fq.ids, pattern = paste0(output_suffix, "_|_demultiplexing_pass.fq.gz"), "")
    write.table(fq.list, fq.list.path, row.names = FALSE, quote = FALSE, col.names = FALSE)
  }

  # 2. Map the reads
  if(premapping == TRUE){
    message("Mapping Ongoing - this will take time")

    minimap2_instances <- ifelse(floor(max_memory/12) < n_cores, floor(max_memory/12), n_cores)
    cl <- parallel::makeCluster(minimap2_instances)
    doParallel::registerDoParallel(cl)
    `%dopar%` <- foreach::`%dopar%`

    minimap2.output.folder <- paste0(output_dir, "/", output_suffix, "_FLASHseqONT_output/", output_suffix, "_mapped_reads")

    # In case of downsampling
    if(is.numeric(downsample_n_reads)){
      minimap2.output.folder <-  paste0(minimap2.output.folder, "_", downsample_n_reads, "/")
    } else {
      minimap2.output.folder <-  paste0(minimap2.output.folder, "/")
    }

    dir.create(minimap2.output.folder)

    invisible(foreach::foreach(x=1:length(fq.list), .packages = "FSnanoporeR", .export = "run_minimap2") %dopar%
      run_minimap2(
        FASTQ = fq.list[x],
        OUTPUT_SUFFIX = fq.ids[x],
        MINIMAP2_Binaries = minimap2_binaries,
        MINIMAP2_REF = minimap2_ref,
        MINIMAP2_BED = minimap2_bed,
        ISOQUANT_MODE = TRUE,
        CHIMERIC_ON = FALSE,
        OUTPUT_DIR = minimap2.output.folder,
        THREADS = ceiling(n_cores/minimap2_instances)
      )
    )
    parallel::stopCluster(cl)

    bam.list <- list.files(minimap2.output.folder, pattern = ".bam$", full.names = TRUE)
    bam.list.path <- paste0(output_dir, "/", output_suffix, "_FLASHseqONT_output/", output_suffix, "_bam.list")
    write.table(bam.list, bam.list.path, row.names = FALSE, quote = FALSE, col.names = FALSE)

  }

  # 4. Run Isoquant
  isoquant.out <- paste0(output_dir, "/", output_suffix, "_FLASHseqONT_output/", output_suffix, "_isoquant")
  # In case of downsampling
  if(is.numeric(downsample_n_reads)){
    isoquant.out <-  paste0(isoquant.out, "_", downsample_n_reads, "")
  }
  if(file.exists(isoquant.out)){unlink(isoquant.out, recursive = TRUE)}

  system(intern = FALSE, command =
           paste0("isoquant.py --read_group file_name --count_exons --clean_start --complete_genedb",
                  " --transcript_quantification ", transcript_counts_type,
                  " --gene_quantification with_inconsistent --data_type nanopore ",
                  " --stranded none",
                  " --check_canonical",
                  " --threads ", n_cores,
                  " --labels ", output_suffix,
                  " --reference ", isoquant_fa,
                  " --genedb ", isoquant_gtf,
                  ifelse(premapping == FALSE, paste0(" --fastq_list ", fq.list.path),  paste0(" --bam_list ", bam.list.path)),
                  " -o ",isoquant.out
           )
  )

}


