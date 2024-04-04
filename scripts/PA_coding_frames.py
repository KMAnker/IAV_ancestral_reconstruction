from Bio import SeqIO

# Open the input file
alignment_file = "ancestral_sequences.fasta"
alignment = SeqIO.parse(alignment_file, "fasta")

# Open the output files
file1 = open("pa.fasta", "w")
file2 = open("pa-x.fasta", "w")

# Define the stop codon
stop_codon = "TAA|TAG|TGA"

# Process each sequence in the alignment
sequences = []
max_length_file1 = 0
max_length_file2 = 0

for seq in alignment:
    seq_id = seq.id
    seq_str = str(seq.seq)
    stop1 = len(seq_str)  # Initialize stop position for file1
    stop2 = len(seq_str)  # Initialize stop position for file2

    # Find the first stop codon in the first reading frame
    for i in range(0, len(seq_str), 3):
        codon = seq_str[i:i+3].upper()
        if codon in stop_codon:
            stop1 = i
            break

    # Find the first stop codon in the new reading frame starting from position 570
    for i in range(571, len(seq_str), 3):
        codon = seq_str[i:i+3].upper()
        if codon in stop_codon:
            stop2 = i
            break

    # Update the maximum length for file1
    max_length_file1 = max(max_length_file1, stop1)

    # Update the maximum length for file2
    max_length_file2 = max(max_length_file2, stop2)

    # Store the sequence and stop positions for file1 and file2
    sequences.append((seq_id, seq_str[:stop1], seq_str[:570] + seq_str[571:stop2]))

# Write to file1
for seq_id, seq1, _ in sequences:
    seq1 += "-" * (max_length_file1 - len(seq1))
    file1.write(">{}\n{}\n".format(seq_id, seq1))

# Write to file2
for seq_id, _, seq2 in sequences:
    seq2 += "-" * (max_length_file2 - len(seq2))
    file2.write(">{}\n{}\n".format(seq_id, seq2))


# Close the output files
file1.close()
file2.close()
