import os
import glob
import anclib

# Define the base directory path
base_path = "PATH to Treetime output directory"
output_base_path = "PATH to output directory"

# List of segments
segments = ["pb2", "pb1", "pa", "ha", "np", "na", "mp", "ns"]

# Loop through each segment folder
for segment in segments:
    segment_path = os.path.join(base_path, segment)
    output_segment_path = os.path.join(output_base_path, segment)
    
    # Create the output segment directory if it doesn't exist
    os.makedirs(output_segment_path, exist_ok=True)
    
    for subdir_timetree in glob.glob(os.path.join(segment_path, "*_timetree")):
        # Extract the segment name from the subdir_timetree
        segment_name = os.path.basename(subdir_timetree).replace("_timetree", "")
    
        # Construct the corresponding mugration subdir path
        subdir_mugration = os.path.join(segment_path, f"{segment_name}_mugration")
    
        # Check if both subdirs exist
        if os.path.exists(subdir_mugration) and os.path.exists(subdir_timetree):
            # Find all .fasta files in the timetree subdir
            fasta_files = glob.glob(os.path.join(subdir_timetree, "*.fasta"))
        
            for fasta_file in fasta_files:
                # Exclude 'ancestral_sequences.fasta'
                if os.path.basename(fasta_file) == "ancestral_sequences.fasta":
                    continue
            
                # Construct file paths
                seq_treefile = os.path.join(subdir_timetree, "timetree.nexus")
                trait_treefile = os.path.join(subdir_mugration, "annotated_tree.nexus")
                trait_probfile = os.path.join(subdir_mugration, "confidence.csv")
            
                # Initialize TreeTimeSeq and TreeTimeTrait objects
                ancseq = anclib.TreeTimeSeq(ancseqfile=fasta_file, seq_treefile=seq_treefile)
                anctrait = anclib.TreeTimeTrait(trait_treefile=trait_treefile, trait_probfile=trait_probfile)
                # Construct AncRecon object
                ancrecon = anclib.AncRecon(ancseq, anctrait)
            
                # Define output file names based on the input file name (without .fasta extension)
                input_filename = os.path.basename(fasta_file)
                input_filename_without_extension = os.path.splitext(input_filename)[0]
                branchdiff_outfile = os.path.join(output_segment_path, f"{input_filename_without_extension}_branchdiffinfo_aa.txt")
                branchdiff_outfile2 = os.path.join(output_segment_path, f"{input_filename_without_extension}_branchdiffinfo_aa_br2.txt")
                branchdiffinfo_outfile3 = os.path.join(output_segment_path, f"{input_filename_without_extension}_branchdiffinfo_aa_br2_all.txt")
                branchinfo_outfile = os.path.join(output_segment_path, f"{input_filename_without_extension}_branchinfo_aa.txt")
                tree_outfile = os.path.join(output_segment_path, f"{input_filename_without_extension}_annotated_tree.nexus")

            
                # Perform the analysis and write output files
                ancrecon.branchdiff(outfilename=branchdiff_outfile, varseq=True, zeroindex=False, printif_traitdiff=False, probmin=0.6, translate=True, br=1)
                ancrecon.branchdiff(outfilename=branchdiff_outfile2, varseq=True, zeroindex=False, printif_traitdiff=False, probmin=0.6, translate=True, br=2)
                ancrecon.branchdiff(outfilename=branchdiffinfo_outfile3, varseq=False, zeroindex=False, printif_traitdiff=False, probmin=0.6, translate=True, br=2)
                ancrecon.branchinfo(outfilename=branchinfo_outfile, varseq=False, zeroindex=False, printif_traitdiff=False, printif_seqdiff=False, probmin=0.6, translate=True, br=1)
                ancrecon.write_tree(outfilename=tree_outfile, label="branchtype", br=1, label_delimiters="braces")
