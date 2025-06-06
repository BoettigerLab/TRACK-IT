import os
import re

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
    junctions = set()  # Use a set to automatically handle duplicates

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
                    junctions.add((ref_name, genome_coord+1, genome_coord + 2)) #moved one base pair downwards because of last CTTCAACTG"T"

                # Search for the reverse complement junction sequence
                rc_junction_pos = read_seq.find(rc_junction_seq)
                if rc_junction_pos != -1:
                    m_number = parseCIGAR(ref_CIGAR) #Aligning to the rightmost end for reverse complement, here I'm using CIGAR string data to do it
                    if m_number>30:
                        genome_coord = ref_start + m_number
                        junctions.add((ref_name, genome_coord-1, genome_coord)) #moved one base pair upwards because of last CTTCAACTG"T"

    except OSError as e:
        print(f"Error reading SAM file: {e}")
        return

    # Sort junctions
    sorted_junctions = sorted(junctions, key=lambda x: (x[0], x[1]))

    # Write junctions to a BED file
    try:
        with open(output_bed, 'w') as bedfile:
            for chrom, start, end in sorted_junctions:
                bedfile.write(f"{chrom}\t{start}\t{end}\n")
    except OSError as e:
        print(f"Error writing BED file: {e}")

# Define the ITR junction sequences
junction_seq = "CTTCAACTGT"
rc_junction_seq = reverse_complement(junction_seq)

# Define the input SAM file and output BED file
samfile_path = r"C:\test\rnd_SB3ITR.sam"
output_bed = r"C:\test\rnd_SB3ITR.bed"

# Find and output the junctions
find_junctions(samfile_path, junction_seq, rc_junction_seq, output_bed)
