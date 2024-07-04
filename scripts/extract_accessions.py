import re
import glob
import sys
import os

# This regular expression matches the "EPI_ISL_number" pattern until the second "|"
gisaid_pattern = re.compile(r'EPI_ISL_[^|]*\|[^|]*')

# Function to extract accessions from a file and add them to two sets
def extract_accessions(file_path):
    # Sets to store unique "EPI_ISL" and other accessions for this file
    gisaid_accessions = set()
    genbank_accessions = set()

    with open(file_path, 'r') as file:
        for line in file:
            columns = line.split('\t')  # split line into columns
            if len(columns) > 1:  # ensure there is a second column
                gisaid_matches = gisaid_pattern.findall(columns[1])  # apply the EPI_ISL pattern to the second column
                if gisaid_matches:
                    gisaid_accessions.update(gisaid_matches)
                else:
                    # If no EPI_ISL pattern is found, take everything before the first underscore
                    accession = columns[1].split('_')[0]
                    if accession:
                        genbank_accessions.add(accession)

    return gisaid_accessions, genbank_accessions

# Check if directory path is provided as command-line argument
if len(sys.argv) < 2:
    print("Please provide the path to the input directory.")
    sys.exit(1)

# Directory to search for input files
directory = sys.argv[1]

# Find all "_namedata.txt" files in the provided directory
input_files = glob.glob(directory + '/*_namedata.txt')

# Extract from all files
for input_file in input_files:
    gisaid_accessions, genbank_accessions = extract_accessions(input_file)

    # Write the unique accessions to a separate output file for each input file
    output_file_path = os.path.basename(input_file).replace('_namedata.txt', '_accessions.txt')
    with open(output_file_path, 'w') as file:
        for accession in gisaid_accessions:
            file.write(accession + '\n')
        for accession in genbank_accessions:
            file.write(accession + '\n')

    print(f'Extracted {len(gisaid_accessions) + len(genbank_accessions)} unique accessions to {output_file_path}')