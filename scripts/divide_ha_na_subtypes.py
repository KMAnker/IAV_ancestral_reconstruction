from Bio import SeqIO

def filter_and_write(input_file, output_file, conditions):
    sequences = []
    with open(input_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            for condition in conditions:
                if condition in record.id:
                    sequences.append(record)
                    break  # Stop checking conditions once a match is found
    SeqIO.write(sequences, output_file, "fasta")

# Process "na.fna" for "n1.fna" and "n2.fna"
filter_and_write("na.fna", "n1.fna", ["H1N1", "H3N1"])
filter_and_write("na.fna", "n2.fna", ["H1N2", "H3N2"])

# Process "ha.fna" for "h1.fna" and "h3.fna"
filter_and_write("ha.fna", "h1.fna", ["H1N1", "H1N2"])
filter_and_write("ha.fna", "h3.fna", ["H3N1", "H3N2"])
