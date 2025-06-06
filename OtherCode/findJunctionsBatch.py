import os
import re
from collections import defaultdict

def parseCIGAR(cigar_string):
    # Use regex to find the integer before 'M'
    match = re.search(r'(\d+)M', cigar_string)
    if match:
        return int(match.group(1))
    else:
        return -1  # or raise an exception, or return a default value

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in reversed(seq))

def find_junctions(samfile_path, junction_seq, rc_junction_seq, output_bed):
    # Use a dictionary to store junctions, their counts, and strandedness
    junction_data = defaultdict(lambda: {'count': 0, 'strand': ''})

    try:
        with open(samfile_path, 'r') as samfile:
            for line in samfile:
                if line.startswith('@'):
                    continue  # Skip header lines

                fields = line.split('\t')
                if len(fields) < 11:
                    continue  # Skip lines that don't have enough fields

                read_seq = fields[9]
                ref_name = fields[2]
                ref_start = int(fields[3]) - 1  # SAM is 1-based, convert to 0-based
                ref_CIGAR = str(fields[5])

                # Search for the forward junction sequence
                junction_pos = read_seq.find(junction_seq)
                if junction_pos != -1:
                    genome_coord = ref_start
                    junction = (ref_name, genome_coord + 1, genome_coord + 2) # Adjust for the sequence offset
                    junction_data[junction]['count'] += 1
                    junction_data[junction]['strand'] = '+'  # Forward strand

                # Search for the reverse complement junction sequence
                rc_junction_pos = read_seq.find(rc_junction_seq)
                if rc_junction_pos != -1:
                    m_number = parseCIGAR(ref_CIGAR)  # Align to the rightmost end using CIGAR string data
                    if m_number > 30:
                        genome_coord = ref_start + m_number
                        junction = (ref_name, genome_coord - 1, genome_coord) # Adjust for the sequence offset
                        junction_data[junction]['count'] += 1
                        junction_data[junction]['strand'] = '-'  # Reverse complement strand

    except OSError as e:
        print(f"Error reading SAM file: {e}")
        return

    # Sort junctions
    sorted_junctions = sorted(junction_data.items(), key=lambda x: (x[0][0], x[0][1]))

    # Write junctions to a BED file with counts and strandedness
    try:
        with open(output_bed, 'w') as bedfile:
            for (chrom, start, end), data in sorted_junctions:
                count = data['count']
                strand = data['strand']
                bedfile.write(f"{chrom}\t{start}\t{end}\t{count}\t{strand}\n")
    except OSError as e:
        print(f"Error writing BED file: {e}")

# Define the ITR junction sequences
junction_seq = "CTTCAACTGT"
rc_junction_seq = reverse_complement(junction_seq)

#Instead uses a folder and looks for the sam files and output junctions accordingly
folderPath = r"C:\test\singleClones"
for file_name in os.listdir(folderPath):
# Check if the file has a .sam extension
    if file_name.endswith('.sam'):
        file_path = os.path.join(folderPath, file_name)
        output_bed = re.sub('.sam', '.bed', file_path)
        find_junctions(file_path, junction_seq, rc_junction_seq, output_bed)
