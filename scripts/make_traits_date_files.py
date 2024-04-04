import sys
import os  # Import the os module to use os.remove()

def extract_traits_and_rename(input_file, output_csv):
    temp_traits_file = "temp_traits.txt"

    with open(input_file, 'r') as infile, open(temp_traits_file, 'w') as temp_outfile:
        for line in infile:
            columns = line.strip().split("\t")
            if any(term in line.lower() for term in ["/swine/", "/sw/", "/pig/", "/sus_scrofa/", "/wild_boar/"]):
                temp_outfile.write(f"{columns[1]}\tswine\n")
            else:
                temp_outfile.write(f"{columns[1]}\thuman\n")

    with open(input_file, "r") as oldnewf, open(temp_traits_file, "r") as tempf, open(output_csv, "w") as outfile:
        outfile.write("name,trait\n")
        renamedict = {words.split()[1]: words.split()[0] for words in oldnewf}
        for line in tempf:
            line = line.rstrip()
            if line:
                for old, new in renamedict.items():
                    if old in line:
                        line = line.replace(old, new)
                outfile.write(line.replace('\t', ',') + "\n")

    os.remove(temp_traits_file)  # Remove the temporary file after it's no longer needed

def extract_dates(input_file, output_csv):
    with open(input_file, 'r') as infile, open(output_csv, 'w') as outfile:
        outfile.write("name,date\n")
        for line in infile:
            parts = line.strip().split("\t")
            second_column = parts[1] if len(parts) > 1 else ""
            if '|' in second_column:
                date = second_column.split('|')[3]
            else:
                split_line = second_column.rsplit('_', 2)
                date = split_line[1] if len(split_line) == 3 else ""
            date = date.rstrip('/').replace('/', '-')
            outfile.write(f"{parts[0]},{date}\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: script.py segment_namedata.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    segment = input_file.split('_namedata.txt')[0]
    traits_csv_path = f"{segment}_traits.csv"
    dates_csv_path = f"{segment}_dates.csv"
    
    extract_traits_and_rename(input_file, traits_csv_path)
    extract_dates(input_file, dates_csv_path)

    print(f"Processing complete. Traits written to {traits_csv_path} and dates to {dates_csv_path}.")
