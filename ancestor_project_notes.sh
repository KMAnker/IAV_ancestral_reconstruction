
################################################################

        ##### Pipeline and notes for the analysis #####

################################################################

### Steps in the pipeline ###

# 1. Collect nucleotide data from NCBI Genbank and GISAID 
    # Downloaded between 26. January and 2. Februrary 2023
    # Only Inf A, HA subtype 1,3; NA subtype 1,2; human and swine sequences downloaded seperately
    # Lists of identifiers for downloaded sequences from GISAID and NCBI Genbank can be found in "data" directory
    # Make sure no spaces are present in Genbank fasta headers: 
        sed -i '' 's/ /_/g' ncbi_human.fna
        sed -i '' 's/ /_/g' ncbi_swine.fna

# 2. Create one file per segment per host:
    divide_by_segment.py

# 3. Concatenate datasets from GISAID and NCBI per segment (still human/swine separately)

# 4. Clean and deduplicate fasta files for each segment
    sequence_cleaner.py
        #settings:
        #  PB2: 2277, 2500, 0
        #  PB1: 2271, 2500, 0
        #  PA: 2148, 2500, 0
        #  HA: 1698, 1900, 0
        #  NP: 1494, 1788, 0
        #  NA: 1407, 1670, 0
        #  MP: 979, 1200, 0
        #  NS: 835, 1100, 0

# 5. Concatenate human and swine datsets per segment 
    # 5.1 Divide HA and NA files per subtype
        divide_ha_na_subtypes.py

# 6. Run cdhit with threshold 97% on all files
    cd-hit-est -i segment.fna -o segment_970.fna -c 0.970 -n 11 -d 0 -M 16000 -T 4

# 7. Re-extract sequences from mixed human-swine clusters
    extract_clusters_cdhit_script.py
    # fasta files were curated using the FLU annotation tool from NCBI (https://www.ncbi.nlm.nih.gov/genomes/FLU/annotation/)
    # sequences with disturbed reading frames, missing start codon, >5 missing bases in the end or mislabeled sequences should be removed

# 8. Align sequences with mafft (auto settings) for each segment (the filtered_segment.fasta files) and inspection in Aliwiew.
   # If a sequence causes large gaps in the rest of the alignment, inspect and maybe remove that sequence
   # Realign and trim the alignment so colums with >40% gaps are removed (https://github.com/agormp/seqconverter) 
        seqconverter -I fasta -O fasta --remfracgapcols 0.4 segment_aln.fasta > segment_aln_trim.fasta

# 9. Rename headers to shorter, unique numbers (https://github.com/agormp/seqconverter)
    seqconverter -I fasta -O fasta --renamenumber seq --savenames namedata.txt segment_aln_trim.fasta > segment_aln_renamed.fasta

# 10. Build tree with IQtree and ancestral reconstuction
    iqtree -s segment_aln_renamed.fasta -pre segment.iqtree -nt auto ntmax 16 -nm 3000 -bb 1000 -asr
    # Output from IQtree can be found in "treefiles>IQtree"

# 11. Make trait and date files ready for treetime
    make_traits_date_files.py
    # 11.1 The treefiles from Iqtree was saved in newick format using Figtree
    # newick-formatted trees, date and traits files can be found in "treefiles>Treetime>inputfiles" 
    # (fasta alignment files are not uploaded)

# 12. We infered a clock tree using treetime while inferring ancestral sequences too
    treetime --tree tree.nwk --dates segment_dates.csv --aln segmetn_aln.fasta --outdir segment_timetree
    # For some segments, using the --covariation option gave a much better tree and model of evolution (PB2, PA)
    # The timetrees and root-to-tip plots were inspected for outliers to see if any sequences should be removed
    # Output from Treetime (except ancestral_sequences.fasta) can be found in "treefiles>Treetime>outputfiles"

# 13. Then, a host annotated tree was inferred with treetime mugration:
    treetime mugration --tree timetree.nexus --states segment_traits.csv --attribute host --confidence --outdir
    # Output from Treetime can be found in "treefiles>Treetime>outputfiles"

# 14. A fasta file was created by Treetime containing the sequences for both tips and ancestral nodes.
    # This file (ancestral_sequences.fasta) was divided into parts coding for each of the viral proteins using
    # custom scripts: 
    "segment"_coding_frames.py

# 15. Using the output files from Treetime, the anclib library was used to summarize sequence and host information
    # in various output files
    anclib_processing.py

    # 15.1. A file containing only the lines for mutated positions in each branch was filtered using
    filter_branchdiff.ipynb

# 16. From the alignment files of the parts coding for specific proteins, amino acid sequence lenghts was 
    # calculated using this script:
    Get_seq_lenghts.py

# 17. Branches and individual positions for each phylogenetic tree and protein-coding alignment were 
    # examined for positive selection pressures using the aBSREL and MEME models from the HyPhy package.
    # The original alignment file used as input for Treetime was converted into seperate files for each
    # protein-coding region and used as input for the hyphy analyses
    "segment"_coding_frames.py
    hyphy absrel --alignment "protein".fasta --tree "segment"_annotated_tree.nexus --output "protein"_absrel.json
    hyphy meme --alignment "protein".fasta --tree "segment"_annotated_tree.nexus --branches "human-swine" --output "protein"_meme_hu-sw.json

        # for the MEME analyses, an annotated treefile with trait categories for each branch (from the anclib output)
        # was used as input, and in four seperate MEME analyses, each of the categories were used as foreground
        # by specifying the "--branches" option

# 18. For the XgBoost analysis, output files from anclib were first formatted and transformed under a one-hot encoding scheme
    file_prep_one_hot_encode.ipynb

    # Then, XgBoost models were trained and evaluated for each protein
    "segment"_xgboost.ipynb


# Most plots and figures were created with R and the respective R code are found in the "figures_tables" directory