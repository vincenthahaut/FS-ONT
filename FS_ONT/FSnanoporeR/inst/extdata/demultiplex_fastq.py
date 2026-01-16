from Bio import SeqIO
import gzip
import csv
import os
#import shutil 
import argparse
from collections import defaultdict

def create_output_folder(outpath='', outsuffix=''):
  output_folder = os.path.join(outpath, f"{outsuffix}_demultiplexing_pass")
  if not os.path.exists(output_folder):
    os.makedirs(output_folder)
  
# Load the read_id and associated umi
def load_umi_orientation_dict(info_path, umi_type='MONOMER', strand_orientated=True):
  orientation_umi_dict = {}
    
  if umi_type == 'MONOMER':
      
    with open(info_path, 'r') as file:
        next(file)  # Skip header line
        for line in file:
            read_id, umi, strand = line.strip().split('\t')
            orientation_umi_dict[read_id] = {'UMI': umi, 'strand': strand}
              
  elif umi_type == 'TRIMER':
      
    with open(info_path, 'r') as file:
        next(file)  # Skip header line
        for line in file:
            read_id,	UMI,	strand,	confidence,	umi_collapsed = line.strip().split('\t')
            orientation_umi_dict[read_id] = {'UMI': umi_collapsed, 'strand': strand}
    
  elif umi_type == 'NONE' and strand_orientated == True:
      
    with open(info_path, 'r') as file:
        next(file)  # Skip header line
        for line in file:
            read_id, barcode, orientation = line.strip().split('\t')
            orientation_umi_dict[read_id] = {'UMI': '', 'strand': orientation}
  
  return orientation_umi_dict



# Load a vsearch object(s)
# Requires: csv
def extract_vsearch_trimming(file_paths=[]):
  data_dict = {}
  for files in file_paths:
      # print(files)
      with open(files, 'r') as file:
          reader = csv.DictReader(file, delimiter='\t')
          #next(reader)
          for row in reader:
            read_id = row['read_id']
            read_section = row['read_section']
            # if the barcode/UMI is at the 5'start of the read
            if read_section == "START":
              bc_end = int(row['read_end'])
              if read_id in data_dict:
                data_dict[read_id]['bc_start_trim'] = max(data_dict[read_id]['bc_start_trim'], bc_end)
              else:
                data_dict[read_id] = {'bc_start_trim': bc_end, 'bc_end_trim': float('-inf')}
            # If the barcode/UMI is at the 3'end of the read
            elif read_section == "END":
              query_l = int(row['query_l'])
              bc_start = int(row['read_start'])
              if read_id in data_dict:
                data_dict[read_id]['bc_end_trim'] = min(data_dict[read_id]['bc_end_trim'], query_l - bc_start)
              else:
                data_dict[read_id] = {'bc_start_trim': float('inf'), 'bc_end_trim': query_l - bc_start }
                    
  return data_dict
  
  
# Index dictionnary
def load_index_unique_dict(detected_bc_path):
  index_unique_dict = {}
  with open(detected_bc_path, 'r') as file:
    reader = csv.DictReader(file, delimiter='\t')
    #next(reader)
    for row in reader:
      read_id = row['read_id']
      index_unique = row['index_unique']
      index_unique_dict[read_id] = index_unique
  return index_unique_dict


# Banned read ids
def load_banlist(file_path):
  banlist_read_ids = []
  with open(file_path, 'r') as file:
    next(file) # skip header
    for line in file:
      read_id = line.strip().split('\t')[0]
      banlist_read_ids.append(read_id)
  return banlist_read_ids

def demultiplex_fastq(fastq_files='', umi_vsearch='', umi_detected='', umi_type='MONOMER', trimming=True, strand_orientated=True, detected_bc='', bc_vsearch='', outpath='', outsuffix='', banlist_path=''):
    # print(strand_orientated)
    # print(trimming)
    # Create a dict to count the number of read associated to each BC
    barcode_counts = defaultdict(int)
    # Create a dict to report information about the reads (length, average phred etc)
    read_length_report = defaultdict(int)
    # Create output folder
    create_output_folder(outpath, outsuffix)
    # Load the umi information and/or strand orientation as dict
    if strand_orientated == 'True' or umi_type in ["MONOMER", "TRIMER"]:
        orientation_umi_dict=load_umi_orientation_dict(umi_detected, strand_orientated=strand_orientated)
    # Load the trimming information as dict.
    if trimming == 'True':
        # print(trimming)
        trim_dict=extract_vsearch_trimming([umi_vsearch, bc_vsearch])
    # Load the banned read ids
    if os.path.isfile(banlist_path):
        #print("Load banlist")
        banlist_read_ids=load_banlist(banlist_path)
    else:
        banlist_read_ids=''
    # Load barcode dictionnary
    indexes_dict=load_index_unique_dict(detected_bc)
    # Create output files for each barcode
    output_files = {}
    for barcode_id in set(indexes_dict.values()):
      outpath_bc=os.path.join(outpath, f"{outsuffix}_demultiplexing_pass")
      output_file_path = os.path.join(outpath_bc, f"{outsuffix}_demultiplexed_{barcode_id}.fq.gz")
      # If previous ones, remove them
      try:
          os.remove(output_file_path)
      except OSError:
          pass
      output_files[barcode_id] = gzip.open(output_file_path, 'wt')
    # Prepare the read length report
    with open(os.path.join(outpath_bc, f"{outsuffix}_read_length.report.txt"), "w") as read_length_report:
        read_length_report.write("read_id\tread_width\tavg_phred_quality\tbarcode_id\tID\n")
        # Read FastQ file record by record
        for fastq_in in fastq_files:
          with gzip.open(fastq_in, 'rt') as fastq_file:
            # limited_fastq_file = islice(fastq_file, 1000)
            # For each reads
            for record in SeqIO.parse(fastq_file, 'fastq'):
                # Extract read ID from the record
                read_id = record.id
                # print(read_id)
                # Exclude reads from banlist and only keep those with a defined barcode
                if read_id not in banlist_read_ids and read_id in indexes_dict:
                    # In case trimming is allowed
                    if trimming == 'True':
                        #print("trimming ongoing")
                        read_length = len(record.seq)
                        if read_id in trim_dict:
                            trim_start = trim_dict[read_id]["bc_start_trim"]
                            trim_end = trim_dict[read_id]["bc_end_trim"]
                            # Adapt trim length to match read_end
                            if trim_end == float('-inf'):
                                trim_end=read_length
                            else:
                                trim_end=read_length-trim_end
                            if trim_start == float('inf'):
                                trim_start=0
                            # Actual trimming
                            record_trimmed = record[trim_start:trim_end]
                            #print('trim-diff:',len(record)-len(record_trimmed))
                        else:
                            record_trimmed = record
                    else:
                      record_trimmed = record
                    # Add UMI and read orientation
                    if umi_type in ['MONOMER', 'TRIMER']:
                        if read_id in orientation_umi_dict:
                            umi_seq=orientation_umi_dict[read_id]['UMI']
                            orientation=orientation_umi_dict[read_id]['strand']
                            record_trimmed.id = f"{read_id}_{umi_seq}"
                            if orientation == '+':
                              # print('test')
                              record_trimmed.seq = record_trimmed.seq.reverse_complement()
                        else:
                            record_trimmed.id = f"{read_id}_NA"
                    # print(record_trimmed.id)
                    # When no UMI is present and orientation is still required
                    if umi_type == 'NONE' and strand_orientated == 'True':
                        #print("strand-oriented")
                        if read_id in orientation_umi_dict:
                            orientation=orientation_umi_dict[read_id]['strand']
                            if orientation == '+':
                                record_trimmed.seq = record_trimmed.seq.reverse_complement()
                    # Demultiplexing
                    barcode_id = indexes_dict[read_id]
                    output_files[barcode_id].write(record_trimmed.format("fastq"))
                    # Increment the number of reads associated to a BC
                    barcode_counts[barcode_id] += 1
                    # Read length report
                    phred_score = sum(record_trimmed.letter_annotations["phred_quality"])
                    read_length = len(record_trimmed.seq)
                    read_length_report.write(f"{record_trimmed.id}\t{read_length}\t{phred_score}\t{barcode_id}\t{outsuffix}\n")

    # Close all output files
    for output_file in output_files.values():
        output_file.close()
    # Save the demultiplexing report
    with open(os.path.join(outpath_bc, f"{outsuffix}_demultiplex.report.txt"), "w") as file:
        file.write("Barcode\tn_reads_before_demultiplexing\tID\n")
        for barcode_id, count in barcode_counts.items():
            file.write(f"{barcode_id}\t{count}\t{outsuffix}\n")

# fastq_in="/home/vincent.hahaut/data_storage/ONT/140224_PAU62267_FSnanopore_760_765/765_FS_FLASHseqONT_output/765_FS.non_chimeric_reads.fq.gz"
# umi_vsearch="/home/vincent.hahaut/data_storage/ONT/140224_PAU62267_FSnanopore_760_765/765_FS_FLASHseqONT_output/765_FS_vsearch.umi.txt"
# umi_type="MONOMER"
# umi_detected="/home/vincent.hahaut/data_storage/ONT/140224_PAU62267_FSnanopore_760_765/765_FS_FLASHseqONT_output/765_FS_detected_umi.txt"
# trimming=True
# strand_orientated=True
# detected_bc='/home/vincent.hahaut/data_storage/ONT/140224_PAU62267_FSnanopore_760_765/765_FS_FLASHseqONT_output/765_FS_detected_barcodes.txt'
# bc_vsearch='/home/vincent.hahaut/data_storage/ONT/140224_PAU62267_FSnanopore_760_765/765_FS_FLASHseqONT_output/765_FS_vsearch.barcodes.txt'
# outpath='/home/vincent.hahaut/data_storage/ONT/140224_PAU62267_FSnanopore_760_765/765_FS_FLASHseqONT_output/'
# outsuffix='765_FS'
# banlist_path="/home/vincent.hahaut/data_storage/ONT/140224_PAU62267_FSnanopore_760_765/765_FS_FLASHseqONT_output/765_FS_chimeric_banList.txt"

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Split a fastq based on a segmentation file')
  parser.add_argument('--input_fastqs_non_chimerics', help='non_chimeric_reads (full path)', required=True)
  parser.add_argument('--input_fastqs_chimerics', help='.chimeric_segments.fq.gz (full path)', required=True)
  parser.add_argument('--umi_vsearch', help='*_vsearch.umi.txt (full path)', required=True)
  parser.add_argument('--umi_detected', help='*_detected_umi.txt (full path)', required=True)
  parser.add_argument('--umi_type', help='MONOMER, TRIMER or NONE', required=True)
  parser.add_argument('--trimming', help='boolean', required=True)
  parser.add_argument('--strand_orientated', help='boolean', required=True)
  parser.add_argument('--bc_vsearch', help='*_vsearch.barcodes.txt (full path)', required=True)
  parser.add_argument('--detected_bc', help='*_detected_barcodes.txt (full path)', required=True)
  parser.add_argument('--outpath', help='Ouput path, typically outsuffix+"_FLASHseqONT_output"', required=True)
  parser.add_argument('--outsuffix', help='Output Suffix', required=True)
  parser.add_argument('--banlist_path', help='*_chimeric_banList.txt (full path)', required=True)  

  args = parser.parse_args()
  
  # Convert to list
  input_fastqs_list = [args.input_fastqs_non_chimerics, args.input_fastqs_chimerics]

  # print(input_fastqs_list)
  demultiplex_fastq(
    fastq_files=input_fastqs_list,
    umi_vsearch=args.umi_vsearch,
    umi_detected=args.umi_detected,
    umi_type=args.umi_type,
    trimming=args.trimming,
    strand_orientated=args.strand_orientated,
    detected_bc=args.detected_bc,
    bc_vsearch=args.umi_vsearch,
    outpath=args.outpath,
    outsuffix=args.outsuffix,
    banlist_path=args.banlist_path)
