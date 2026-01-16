import pysam
import argparse
import concurrent.futures 
from multiprocessing import Manager, Pool

# Function to call add_tags_to_bam
def process_bam(args):
    bam_in, tag_data, cell_barcode, bam_out = args
    add_tags_to_bam(bam_in, bam_out, tag_data, cell_barcode)

# Using a dict with read_id as key and new tags as items, update the tags of a BAM file
def add_tags_to_bam(bam_in, bam_out, tag_data, cell_barcode):
    with pysam.AlignmentFile(bam_in, "rb") as bam, pysam.AlignmentFile(bam_out, "wb", template=bam) as output_bam:
            for read in bam:
                read_id = read.query_name
                if read_id in tag_data:
                    tags = tag_data[read_id]
                    for tag_id, tag_value in tags.items():
                        read.set_tag(tag_id, tag_value)
                read.set_tag("CB", cell_barcode)
                output_bam.write(read)
    pysam.index(bam_out)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add new SAM tags to BAM files')
    parser.add_argument('--threads', help='Number of BAM files to run in parallel', required=True)
    parser.add_argument('--bam_in', nargs='+', help='list of input bam (full path)', required=True)
    parser.add_argument('--bam_out', nargs='+', help='list of output bam (full path)', required=True)
    parser.add_argument('--cell_barcodes', nargs='+', help='list of CB values associated with BAM files', required=True)
    parser.add_argument('--tag_file', help='Path to the tag file - expects read_id in first column and new tags in separated columns (tag name as column.name)', required=True)

    args = parser.parse_args()

    # Create a Manager to manage shared data (tag_data is quite heavy)
    # Avoid sending tag_data to each worker
    #manager = Manager()
    #tag_data_shared = manager.dict()

    # Read the tag file and store tags for each read ID
    print("1. Load Tag Data")
    print("ThreadPool used")
    #tag_data = {}
    #with open(args.tag_file, 'r') as tags:
    #    header = tags.readline().strip().split()  # Read the header to use as tag IDs
    #    for line in tags:
    #        parts = line.strip().split()
    #        read_id = parts[0]
    #        tags_data = dict(zip(header[1:], parts[1:]))  # Use header as tag IDs and create tag_id:tag_value pairs
    #        tag_data_shared[read_id] = tags_data

    tag_data = {}
    with open(args.tag_file, 'r') as tags:
        header = tags.readline().strip().split()  # Read the header to use as tag IDs
        for line in tags:
            parts = line.strip().split()
            read_id = parts[0]
            tags_data = dict(zip(header[1:], parts[1:]))  # Use header as tag IDs and create tag_id:tag_value pairs
            tag_data[read_id] = tags_data

    # Run in parallel
    #bam_args = [(bam, tag_data_shared, cell_barcode, bam_out) for bam, cell_barcode, bam_out in
    #            zip(args.bam_in, args.cell_barcodes, args.bam_out)]
    bam_args = [(bam, tag_data, cell_barcode, bam_out) for bam, cell_barcode, bam_out in
                zip(args.bam_in, args.cell_barcodes, args.bam_out)]

    print("2. Start Processing")
    with concurrent.futures.ThreadPoolExecutor(max_workers=int(args.threads)) as executor:
        results = executor.map(process_bam, bam_args)

    #with Pool(processes=int(args.threads)) as pool:
    #    results = pool.map(process_bam, bam_args)

    # Clean-up in case
    # del tag_data_shared
    print("3. Finished!")



# bam_in=["/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_mapped_reads_50000//BC129_sorted.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_mapped_reads_50000//BC130_sorted.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_mapped_reads_50000//BC131_sorted.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_mapped_reads_50000//BC134_sorted.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_mapped_reads_50000//BC136_sorted.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_mapped_reads_50000//BC138_sorted.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_mapped_reads_50000//BC144_sorted.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_mapped_reads_50000//BC146_sorted.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_mapped_reads_50000//BC147_sorted.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_mapped_reads_50000//BC150_sorted.bam"]
# cell_barcodes=["BC129","BC130","BC131","BC134","BC136","BC138","BC144","BC146","BC147","BC150"]
# tag_file="/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_isoquant_50000/barcode01/barcode01_read_assignments.updated.txt"
# bam_out=["/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_isoquant_50000/barcode01//umi/barcode01_BC129_isoquant_umi.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_isoquant_50000/barcode01//umi/barcode01_BC130_isoquant_umi.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_isoquant_50000/barcode01//umi/barcode01_BC131_isoquant_umi.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_isoquant_50000/barcode01//umi/barcode01_BC134_isoquant_umi.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_isoquant_50000/barcode01//umi/barcode01_BC136_isoquant_umi.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_isoquant_50000/barcode01//umi/barcode01_BC138_isoquant_umi.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_isoquant_50000/barcode01//umi/barcode01_BC144_isoquant_umi.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_isoquant_50000/barcode01//umi/barcode01_BC146_isoquant_umi.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_isoquant_50000/barcode01//umi/barcode01_BC147_isoquant_umi.bam","/home/vincent.hahaut/data_storage/ONT/180823_PAQ21221_FSnanopore//barcode01_FLASHseqONT_output/barcode01_isoquant_50000/barcode01//umi/barcode01_BC150_isoquant_umi.bam"]

