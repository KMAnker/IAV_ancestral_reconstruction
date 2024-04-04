from Bio import SeqIO

# Open the input file
alignment_file = "ancestral_sequences.fasta"
alignment = SeqIO.parse(alignment_file, "fasta")

# Open the output files
file1 = open("n1.fasta", "w")

# Define the stop codon
stop_codon = "TAA|TAG|TGA"

# Find the first stop codon in the first reading frame for each sequence
stops = {}
for seq in alignment:
    for i in range(0, len(seq.seq)-2, 3):
        codon = str(seq.seq[i:i+3]).upper()
        if codon in stop_codon:
            stop = i
            break
    else:
        stop = len(seq.seq)
    stops[seq.id] = stop

# Reset the file pointer and iterate through the input file again
alignment = SeqIO.parse(alignment_file, "fasta")
for seq in alignment:
    # Write the first file
    seq1 = seq.seq[0:stops[seq.id]]
    seq1 += "-"*(max(stops.values())-len(seq1))
    file1.write(">{}\n{}\n".format(seq.id, seq1))

# Close the output file
file1.close()
