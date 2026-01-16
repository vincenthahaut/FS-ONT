import pandas as pd
import argparse
import gzip

def parse_segment_file(segment_file):
    """Parse the segment file and return a DataFrame with selected columns."""
    selected_columns = ["read_id", "segment_start", "segment_end", "putative_chimeric", "fastq_id_new"]
    segment_df = pd.read_csv(segment_file, sep='\t', usecols=selected_columns)
    # Filter rows where putative_chimeric is TRUE
    segment_df = segment_df[segment_df["putative_chimeric"] == True]
    # Transform to dictionary
    read_dict = {}
    for index, row in segment_df.iterrows():
        read_id = row["read_id"]
        segment_start = row["segment_start"]
        segment_end = row["segment_end"]
        fastq_id_new = row["fastq_id_new"]
        #
        if read_id not in read_dict:
            read_dict[read_id] = []
        # Append the new_read_id to the list of new_read_ids for the corresponding read_id
        read_dict[read_id].append({'fastq_id_new': fastq_id_new, 'segment_start': segment_start, 'segment_end': segment_end })
    return read_dict


def split_fastq(input_fastq, segment_dict, output_fastq_chimeric):
    """Split FASTQ by specified positions."""
    with gzip.open(input_fastq, 'rt') as input_fh, gzip.open(output_fastq_chimeric, 'wt') as output_fastq_chimeric_w:
        for i, line in enumerate(input_fh):
            if (i // 4) % 4 == 0:
                # Extract sequence/quality - adjusting for R/python base system
                read_id = line.strip().split(' ')[0][1:]
                sequence = input_fh.readline()
                _ = input_fh.readline()  # Skip the "+" line
                quality = input_fh.readline()
                # If Chimeric
                if read_id in segment_dict:
                    for segment_info in segment_dict[read_id]:
                        new_fastq_id = segment_info['fastq_id_new']
                        segment_start = segment_info['segment_start']
                        segment_end = segment_info['segment_end']
                        # Split reads
                        sequence_split = sequence.strip()[segment_start-1:segment_end-1]
                        quality_split = quality.strip()[segment_start-1:segment_end-1]
                        # Write modified sequence with new read_id
                        output_fastq_chimeric_w.write(f"@{new_fastq_id}\n{sequence_split}\n+\n{quality_split}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Split a fastq based on a segmentation file')
    parser.add_argument('--input_fastq', help='input.fastq file (full path)', required=True)
    parser.add_argument('--segment_file', help='chimericSegments.info.txt (full path)', required=True)
    parser.add_argument('--output_fastq_chimeric', help='output.fastq chimeric segments post-split (full path)', required=True)

    args = parser.parse_args()

    # Parse the segment file
    segment_dict = parse_segment_file(args.segment_file)

    # Process the FASTQ file and write to the output file
    split_fastq(args.input_fastq, segment_dict, args.output_fastq_chimeric)


