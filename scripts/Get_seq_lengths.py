import sys

# Define a function to calculate sequence length without gaps
def calculate_sequence_length(sequence):
    length = len(sequence.replace('-', ''))
    return int(length / 3) if length % 3 == 0 else length / 3

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 3:
    print("Usage: python Get_seq_lengths.py input.fasta output.txt")
    sys.exit(1)

# Input and output file paths from command-line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Initialize variables to store header and sequence
header = ''
sequence = ''
sequence_lengths = []

# Open the input FASTA file for reading
try:
    with open(input_file, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                # If a header line is encountered, save the previous sequence length (if any)
                # and store the new header
                if header:
                    sequence_lengths.append((header, calculate_sequence_length(sequence)))
                header = line[1:]  # Remove the '>' character
                sequence = ''  # Initialize the sequence
            else:
                # Concatenate sequence lines
                sequence += line

        # Calculate and store the length of the last sequence in the file
        if header:
            sequence_lengths.append((header, calculate_sequence_length(sequence)))

    # Write the header and sequence lengths to the output file
    with open(output_file, 'w') as output:
        for header, length in sequence_lengths:
            output.write(f'{header}\t{length}\n')

    print(f'Sequence lengths have been written to {output_file}')
except FileNotFoundError:
    print(f"Error: File '{input_file}' not found.")
    sys.exit(1)
except Exception as e:
    print(f"An error occurred: {e}")
    sys.exit(1)
