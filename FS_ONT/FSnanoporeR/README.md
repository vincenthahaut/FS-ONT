# FLASH-seq-ONT - FSnanoporeR Manual

This repository contains the scripts associated to: [address]

The main part of the code is contained in a R package ('FSnanoporeR') which provides all the necessary functions to process FLASH-seq data sequenced with Oxford Nanopore long-read sequencing technology (ONT).


## Dependencies

Create a conda environment with the required dependencies:

```
mamba env create -n FSnanopore --file conda_environment.yaml
```

We strongly recommend using mamba instead of conda to manage this environement.

## FSnanoporeR package Installation

### Method-1

Download the github R package from: [FSnanopore](https://github.com/vincenthahaut/IOB/tree/master/6_Nanopore-FS/FSnanoporeR)

Active the conda environment:

```
conda activate FSnanopore
```

Run the following command:

```
R CMD INSTALL /path/to/R/package/FSnanoporeR_folder
```

or use the pre-build binary:

```
R CMD INSTALL FSnanoporeR_0.1.0.tar.gz
```

### Method-2 (not tested)

```
conda activate FSnanopore

conda skeleton cran https://github.com/vincenthahaut/IOB/tree/master/6_Nanopore-FS/FSnanoporeR

R --version

conda build --R=4.2.3 r-fsnanoporer
```


## Preparing the input

This pipeline requires basecalled FASTQ files (dorado or guppy) demultiplexed (guppy) by plate index. Follow the latest guidelines from ONT as dorado / guppy tend to be often updated. 

```
# OUTDATED: If you asked MinKnow (older versions) to output fast5 and not pod5
pod5 convert fast5 ./fast5/*.fast5 --output pod5/ --one-to-one ./fast5/

# Basecalling (command line or on minKNOW)
# NB: Use at least high-accuracy models!
dorado_model="dna_r10.4.1_e8.2_400bps_sup@v4.2.0"
dorado basecaller "$dorado_model" pod5/ > unmapped_reads_with_moves.bam

# Convert BAM to fastq
samtools fastq unmapped_reads_with_moves.sam > unmapped_reads_with_moves.fastq

# Split FASTQ into smaller chunks
# Improves the performance of guppy_barcoder which crashes if fed a too big single file.
mkdir fastq
cd fastq
split -l 1000000 ../unmapped_reads_with_moves.fastq 
for file in x*; do mv "$file" "$file".fastq; done

# Demultiplexing by Plate Index
cd ..
THREADS=48
guppy_barcoder 
--fastq_out 
--barcode_kits SQK-NBD114-24 
--min_score_adapter 30 
--min_score_barcode_rear 30 
--min_score_barcode_front 30 
--compress_fastq 
-t $THREADS
-i ./fastq/ 
-s demultiplexed_30adapterscores/
```

## Pipeline

### Rscript Wrapper

Presently, the entire pipeline can be run at once using a R wrapper script.

The wrapper can be found in the R package at 'FSnanoporeR/inst/extdata/wrapper_FS.R' and requires a YAML samplesheet as an argument (example FSnanoporeR/inst/extdata/intialise_FS_ONT_example.yaml").

To run it from the terminal, in a linux environment, use:

```
# Active the conda environment
conda activate FSnanopore

# Run the wrapper
ulimit -n 100000
Rscript FSnanoporeR/inst/extdata/wrapper_FS.R sample_sheet.yaml
```

**Note:** This script relies on several tools (blast, vsearch, seqkit,...). Running the R wrapper from the terminal in an activated conda session will call them from your conda folder. This wrapper cannot be called from Rstudio as this GUI works from /usr/bin/R.


### YAML samplesheet

The YAML samplesheet contains a series of options that can be set. An example of the samplesheet, with their description, can be found at 'FSnanoporeR/inst/extdata/intialise_FS_ONT_example.yaml'.

**Note:** If running Isoquant, please use Gencode annotation or manually update the isoquant option '--complete_genedb' accordingly (hardcoded).

## Speed, Parallel Processing and IT requirements

This pipeline is computationally intensive, especially during the mapping steps, read assignment to a feature and recovery of inconsistent reads. The most RAM hungry R functions are progressively replaced with python streaming versions.

To speed up processing, the provided FASTQs are divided into chunks of equal read numbers. Chunks are processed in parallel on a single core. Both the number of cores and chunks can be defined in the YAML file. Same safeguards have been implemented to cap the number of processes run in parallel and avoid crashes.

A full promethion flowcell will generate folders of processed >400Gb of data. The pipeline should not be called with <800 Gb of storage space. 

Pipeline was developed on Ubuntu 20.04, 204 Gb RAM, 32 cores. On a full dataset (>1.9Tb raw pod5 basecalled), the pipeline takes ~24-48hrs to run depending on the amount of data and percentage of chimeric reads. At this stage, I do not recommend running it on a device with <100 Gb RAM. The next step for us is to optimise the RAM cosumption / storage.


## Pipeline Description

The script can be divided into 14 tasks some which are optional. 

1. Split reads in chunks for parallel processing.
2. Check for presence of an ISPCR inside the read (= chimeric) using BLAST. Split reads into segments if found.
3. **Optional:** Map chimeric reads and segments.
4. **Optional:** Validate segmentation by comparing the alignments of chimeric segments and reads.
5. Extract cell barcodes positions using VSEARCH.
6. Extract and correct cell barcode sequences.
7. Extract the UMI (monomeric or trimeric) or the position of the oligo-dT using VSEARCH.
7. Demultiplex reads:
 - **Optional:** Trim reads to remove the barcodes / PCR adapters.
 - **Optional:** Correct read orientation based on UMI or oligo-dT orientation.
8. Combine files generated during parallel processing.
9. **Optional:** Generate an HTML report of the previous tasks.
10. **Optional:** Cleanup intermediate files.
11. **Optional:** Map reads.
11. **Optional:** Run Isoquant on mapped reads and assign them to a gene/transcript features.
12. **Optional - If UMI**: Convert isoquant output to a UMI_tools compatible BAM file.
13. **Optional - If UMI**: Run UMI_tools to counts the UMI associated to a gene- or transcript-feature.
14. **Optional - If UMI**: Combine UMI_tools results into a matrix of count (gene- or transcript-).

## Outputs

The pipeline generates a series of files. The following paragraph describes the different outputs. Depending on your selected options, some files may be absent.

With OUTPUT_DIR and OUTPUT_SUFFIX as arguments from the samplesheet:

```
└── OUTPUT_FOLDER
    ├── OUTPUT_SUFFIX_demultiplexing_pass                       # Folder containing the demultiplexed FASTQ
    ├── OUTPUT_SUFFIX.chimeric_reads.fq.gz                      # FASTQ file containing the suspected chimeric reads pre-segmentation
    ├── OUTPUT_SUFFIX.chimeric_segments.fq.gz                   # FASTQ file containing the suspected chimeric reads post-segmentation
    ├── OUTPUT_SUFFIX.non_chimeric_reads.fq.gz                  # FASTQ file without suspected chimeric reads
    ├── OUTPUT_SUFFIX_bam.list                                  # List of the BAM files path for isoquant
    ├── OUTPUT_SUFFIX_blast_search.txt                          # Blast search of the PCR adapter sequence
    ├── OUTPUT_SUFFIX_chimericReads.ids.txt                     # Read IDs of the suspected chimeric
    ├── OUTPUT_SUFFIX_chimericSegments.info.txt                 # Summary of the blast search at the read-level. Contains a description of the segment for each read, chimeric or not. 
    ├── OUTPUT_SUFFIX_chimeric_banList.txt                      # Read IDs of the chimeric reads which display problematic features (bad segmentation and/or in-silico)
    ├── OUTPUT_SUFFIX_demultiplex.report.txt                    # Number of reads, median/sd length associated to each cell barcode
    ├── OUTPUT_SUFFIX_detected_barcodes.txt                     # Summary of the read-associated barcodes, post-filtering
    ├── OUTPUT_SUFFIX_detected_barcodes.unfiltered.txt          # Summary of the read-associated barcodes, pre-filtering
    ├── OUTPUT_SUFFIX_detected_umi.txt                          # Read associated UMI
    ├── OUTPUT_SUFFIX_extractBarcodes.report.txt                # Various information about the detected barcodes, per chunk (#, mismatches, perfect, ...)
    ├── OUTPUT_SUFFIX_fq.list                                   # Full path of the FASTQ reads post-demultiplexing - for isoquant
    ├── OUTPUT_SUFFIX_get_UMI.report.txt                        # Various information about the detected UMI, per chunk (#, masked bases, mismatches, ...)
    ├── OUTPUT_SUFFIX_in_silico_chimerics.log.txt               # Summary of the suspected chimeric segments mapping. Comparison of the segments vs segment overlaps.
    ├── OUTPUT_SUFFIX_mapping_chimerics.log.txt                 # Summary of the suspected chimeric segments mapping. Comparison of the segments vs reads overlaps.
    ├── OUTPUT_SUFFIX_chimeric_sorted.bam                       # Minimap2 alignment of the chimeric reads and segments altogether
    ├── OUTPUT_SUFFIX_read_length.report.txt                    # Read ID, length and sum of the base quality.
    ├── OUTPUT_SUFFIX_report.html                               # HMTL report of the demultiplexing
    ├── OUTPUT_SUFFIX_splitChimeric_reads_report.txt            # Vsearch results of the cell barcode
    ├── OUTPUT_SUFFIX_vsearch.barcodes.txt                      # Vsearch results of the cell barcode
    ├── OUTPUT_SUFFIX_vsearch.umi.txt                           # Vsearch results of the UMI detection
    ├── OUTPUT_SUFFIX_mapped_reads                              # Mapped reads (sorted.bam and index)
    └── OUTPUT_SUFFIX_isoquant                                  # Folder containing the results of the isoquant pipeline
        └── OUTPUT_SUFFIX                               
            ├── umi                                                         # Folder containing Isoquant BAM files converted to UMI_tools compatible BAMs, per cell
            ├── umi_tools                                                   # Folder containing the UMI_tools quantification, per cell
            ├── aux                                                         # Empty folder created by isoquant
            ├── OUTPUT_SUFFIX_isoquant_umitools_counts.gene.tsv             # UMI Matrix of counts, gene-level
            └── OUTPUT_SUFFIX_isoquant_umitools_counts.transcripts.tsv.gz   # UMI Matrix of counts, transcript-level

```

## Recommandation for Post-Processing

Post-processing of FLASH-seq-ONT results can be performed with any pipeline dedicated to long-read single-cell RNA-sequencing. If not using the pipeline incorporated in FSnanoporeR package, mapping and read count can be performed as follows:

### Mapping

I recommend using minimap2. However, this mapper is so good that it can even map confidently reads from a different species on the human genome. If you are after high quality mapping, I recommend adding an extra filter on gap-compressed per-base sequence divergence or 'de' flag.

```
# Reference
minimap2 -d genome_reference.mni /path/to/genome_reference.fa.gz
paftools.js gff2bed /path/to/genome_annotation.gtf.gz > genome_annotation.bed

# Standard
minimap2 -ax splice -I12g -ub -t 5 --junc-bed genome_annotation.bed genome_reference.mni $file > aligned_files.bam
samtools view -F 260 -b aligned_files.bam | samtools sort -@ 5 -o aligned_files.sorted.bam
samtools index aligned_files.sorted.bam

# Highest quality
samtools view -F 260 -e "[de] < 0.075" -b aligned_files.bam | samtools sort -@ 5 -o aligned_files.sorted.bam
```

### Read Counting

Both isoquant and Bambu work well on FLASH-seq-ONT data. However, the additional outputs from isoquant on feature mappings and its lower false positive rate (ref[https://www.nature.com/articles/s41587-022-01565-y]) are the reasons why it is included in the main pipeline. 

Some of the options from isoquant are going beyond the scope of the demonstration of FLASH-seq-ONT. Therefore, the R wrapper includes premade choices such as the absence of novel transcript detection, gene counting or transcript counting. To run isoquant with your parameters on the FLASH-seq-ONT results, modify directly the run_isoquant.R function or update the following command with your arguments:

```
FASTQ=OUTPUT_SUFFIX_fq.list     # Full path of the FASTQ to process from the FSnanoporeR pipeline

isoquant.py 
--read_group file_name \
--count_exons --clean_start \ 
--complete_genedb \                             # Only valid for gencode annotation
--transcript_quantification unique_only \       # Could be set differently if you are confident about your reads
--gene_quantification with_inconsistent \ 
--data_type nanopore \ 
--stranded none \                               # This option currently does not influence isoquant read assignment
--threads 20 \
--labels $OUTPUT_SUFFIX \
--reference /path/to/reference.fa \
--genedb /path/to/reference.gtf \
--fastq_list $FASTQ \
-o $OUTPUT_DIR

# --report_novel_unspliced could be added
```


**Notes:** Using different different options in isoquant may results in unforseen effects on the umi_tools BAM conversion. 


### Update Isoquant count matrices

Isoquant is a great tool to count reads associated to certain genes / transcripts. Unfortunately, single-cell RNA-sequencing data present some challenges that pushed us to revise some of the read counting strategies elaborated in Isoquant.

Isoquant defines different levels of certainty regarding the association of a read to a feature. We first select:

```
unique
unique_minor_difference
inconsistent
```

The mapping inconsistencies are also reported for each read. We have observed that some of them are more prone to be associated to problematic mappings:

```
# Introns
"intron_migration"
"intron_alternation"
# Exons
"major_exon_elongation"
"mutually_exclusive_exons"
"exon_skipping"
"exon_merge"
"exon_gain"
"exon_detach"
# Splicing
"alt_donor_site_novel"
"alt_acceptor_site_novel"
# Other
"alternative_structure"
# Typical from highly truncated molecules / gDNA
"mono_exonic"
"antisense"
"incomplete_intron_retention"
"internal_polya"
```

We exclude the reads that display one of these features except for the uniquely assigned reads with the "mono_exonic_match" (= good match to a single-exon gene). This allows us to refine the categories to:

```
unique
unique_minor_difference
inconsistent_recovered
```

Using Isoquant *.read_assignments.tsv file, we update the BAM file with additional tags (see R ```update_read_assignments```):

```
XT - Isoform id - Only consider uniquely assigned reads (unique / unique_minor_difference)
TS - XT associated assignment status - 'Assigned' or 'Unassigned'
XI - Isoform id - Consider uniquely assigned and recovered inconsistent reads 
IS - XI associated assignment status - 'Assigned' or 'Unassigned'
XG - Gene id - Consider all uniquely assigned and inconsistent reads (= isoquant with_inconsistent)
GS - GS associated assignment status - 'Assigned' or 'Unassigned'
XR - Gene id - Consider uniquely assigned and recovered inconsistent reads 
RS - RS associated assignment status - 'Assigned' or 'Unassigned'
CB - cell barcode
UB - UMI
IF - Isoquant Flag(s)
IC - Isoquant classification
IN - Isoquant canonical splicing
```

Finally, we also correct one decision from Isoquant to count multiple times reads assigned to several features. Any read assigned to more than one feature (at the gene- or isoform-level) is discarded.
