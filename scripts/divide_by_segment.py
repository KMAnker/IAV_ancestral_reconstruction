from Bio import SeqIO
import os
import sys

def main():
    if len(sys.argv) != 2:
        print("Usage: python divide_by_segment.py input_file")
        return
    
    input_file = sys.argv[1]
    
    # Iterate over the sequences in the input file
    for seq_record in SeqIO.parse(input_file, "fasta"):
        header = seq_record.id

        if header.endswith("|1") or header.endswith("_1"):
            output_file = input_file.replace(".fna", "_pb2.fna")
        elif header.endswith("|2") or header.endswith("_2"):
            output_file = input_file.replace(".fna", "_pb1.fna")
        elif header.endswith("|3") or header.endswith("_3"):
            output_file = input_file.replace(".fna", "_pa.fna")
        elif header.endswith("|4") or header.endswith("_4"):
            output_file = input_file.replace(".fna", "_ha.fna")
        elif header.endswith("|5") or header.endswith("_5"):
            output_file = input_file.replace(".fna", "_np.fna")
        elif header.endswith("|6") or header.endswith("_6"):
            output_file = input_file.replace(".fna", "_na.fna")
        elif header.endswith("|7") or header.endswith("_7"):
            output_file = input_file.replace(".fna", "_mp.fna")
        elif header.endswith("|8") or header.endswith("_8"):
            output_file = input_file.replace(".fna", "_ns.fna")
        else:
            output_file = input_file.replace(".fna", "_extra.fna")

        # Open the output file for writing
        with open(output_file, "a") as output_handle:
            # Write the sequence to the output file
            SeqIO.write(seq_record, output_handle, "fasta")

if __name__ == "__main__":
    main()
