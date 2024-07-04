# This code comes from Biopython: https://biopython.org/wiki/Sequence_Cleaner
# It is modified in order to have one more parameter that remove sequences that are longer that the specified max_length (3rd input parameter)
# and other features have been added, like checking that sequences only contain allowed characters, that they can max have 5 ambiguous characters, 
# that fasta headers don't have spaces or other characters except those allowed by IQtree, and that sequences from environmental or synthetic samples 
# are removed as well as removal of identical sequences so there are no duplicates.

import sys
import os
import re
from Bio import SeqIO


def sequence_cleaner(fasta_file, min_length=0, max_length=0, por_n=100):
    # Create our hash table to add the sequences
    sequences = {}
    # Change path
    abspath = os.path.abspath(".")
    bn = os.path.basename(fasta_file)
    file_path = os.path.dirname(fasta_file)
    output_path = abspath
    os.chdir(file_path)

    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(bn, "fasta"):
        # Take the current sequence
        sequence = str(seq_record.seq).upper()
        # Check if the current sequence is according to the user parameters
        if (
            len(re.findall("[ACGTNRYKMSWBDHV-]", sequence)) == len(sequence)
            and len(sequence) >= min_length
            and (float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n
            and len(sequence) <= max_length
            and len(re.findall("[RYKMSWBDHV]", sequence)) < 5

        ):
            seq_record.id = re.sub(r'[\s?,:;+\'()]', '_', seq_record.id)
            if re.search("environment", seq_record.id, re.IGNORECASE):
                continue
            if re.search("synthetic", seq_record.id, re.IGNORECASE):
                continue      
            # If the sequence passed in the test "is it clean?" and it isn't in the
            # hash table, the sequence and its id are going to be in the hash
            if sequence not in sequences:
                sequences[sequence] = seq_record.id
            # If it is already in the hash table, we will not include it in the final fasta output file
            else:
                continue

    # Write the clean sequences

    # Create a file in the same directory where you ran this script

    with open(os.path.join(output_path,("clear_" + bn)), "w+") as output_file:
        # Just read the hash table and write on the file as a fasta format
        for sequence in sequences:
            output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")

    print("CLEAN!!!\nPlease check clear_" + bn)


userParameters = sys.argv[1:]


try:
    if len(userParameters) == 1:
        sequence_cleaner(userParameters[0])
    elif len(userParameters) == 2:
        sequence_cleaner(userParameters[0], float(userParameters[1]))
    elif len(userParameters) == 3:
        sequence_cleaner(userParameters[0], float(userParameters[1]), float(userParameters[2]))
    elif len(userParameters) == 4:
        sequence_cleaner(userParameters[0], float(userParameters[1]), float(userParameters[2]), float(userParameters[3]))
    else:
        print("There is a problem!")
except:
    print("There is a problem!")

