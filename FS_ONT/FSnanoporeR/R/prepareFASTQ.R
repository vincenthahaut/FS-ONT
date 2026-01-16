#' @title Split FASTQ file by line and filter by length
#'
#' @description
#' Split a FASTQ file in equivalent chunks of reads based on the number of threads.
#'
#' @param FASTQ character. "/path/to/the/input/fastq.gz". Gziped (expect '.gz') or not.
#' @param OUTPUT_DIR character. "/path/to/output/"
#' @param OUTPUT_SUFFIX character. sample id.
#' @param N_CHUNKS numeric. In how many pieces should the FASTQ be split ?
#' @param THREADS numeric.
#' @param silent boolean. If TRUE will print various messages.
#'
#' @importFrom plyr round_any
#'
#' @export
prepareFASTQ <- function(FASTQ = NULL,
                         OUTPUT_DIR = NULL,
                         N_CHUNKS = NULL,
                         THREADS = NULL,
                         OUTPUT_SUFFIX = NULL,
                         silent = FALSE){
  
  options(scipen=999)
  
  if(file.exists(paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_splitFASTQ/"))){unlink(paste0(OUTPUT_DIR, "/splitFASTQ/"),recursive=TRUE)}
  if(!file.exists(OUTPUT_DIR)){dir.create(OUTPUT_DIR)}
  dir.create(paste0(OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_splitFASTQ/"))
  
  # If provide a single FASTQ file - old approach from guppy
  if(!file.info(FASTQ)$isdir){
    if(silent){message("Provided a single FASTQ - split in smaller files")}
    if(silent){message("This may take some time - better provide the folder resulting from guppy demultiplexing in a --one-to-one fashion")}
    # Calculate ideal number of lines to split
    # LINES / (N_cores / 4)
    # ==> Each subset will be processed with 4 cores
    if(grepl(".gz", FASTQ)){
      LINES <- as.numeric(system(intern = TRUE, command = paste0("zcat ", FASTQ, "| wc -l ")))/4
    } else {
      LINES <- as.numeric(system(intern = TRUE, command = paste0("wc -l ", FASTQ)))/4
    }
    if(!silent){message(paste0("Detected: ", LINES, " reads"))}
    
    LINES_SUBSET <- plyr::round_any((LINES*4)/((N_CHUNKS)-1), as.numeric(paste0(1, paste0(rep(0,nchar(floor(LINES/N_CHUNKS-1)),collapse = ""),collapse = ""))))
    
    if(!silent){message(paste0("Split by ", LINES_SUBSET/4, " lines into ", (N_CHUNKS)-1, " files"))}
    
    # Split file
    system(intern = FALSE,
           command = paste0("zcat ", FASTQ, "| split --additional-suffix=.fq -l ", LINES_SUBSET, " - ", OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_splitFASTQ/", OUTPUT_SUFFIX, "_splitFASTQ")
    )
    
    # Else if providing a directory
  } else if(file.info(FASTQ)$isdir){
    if(!silent){message("Directory provided")}
    
    fastq_list <- list.files(FASTQ, full.names = TRUE, pattern = ".fastq")
    
    num_elements <- length(fastq_list)
    elements_per_vector <- ceiling(num_elements / N_CHUNKS)
    
    fastq_list_split <- split(fastq_list, rep(1:N_CHUNKS, each = elements_per_vector, length.out = num_elements))
    
    cmd <- ifelse(all(grepl(".gz", fastq_list)), "zcat", "cat")
    
    parallel::mclapply(1:length(fastq_list_split), mc.silent = silent, mc.cores = THREADS, function(x)
      
      system(intern = FALSE,
             command = paste0(cmd, " ", paste0(fastq_list_split[[x]], collapse = " "), " > ", OUTPUT_DIR, "/", OUTPUT_SUFFIX, "_splitFASTQ/", OUTPUT_SUFFIX, "_splitFASTQ", x ,".fq")
      )
      
    )
    
    return("Finished")
  }
  
}
