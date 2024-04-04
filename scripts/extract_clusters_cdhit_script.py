
import sys
import os
#specify input file from cd-hit (clstr file)
inputfile = sys.argv[1]
#specify the original (cleaned) fasta file 
inputfile2 = sys.argv[2]
#name of the filtered cd-hit clstr file
outputfile = "filtered_" +  inputfile
#name of the headers in the filtered cd-hit file
outputfile2 = os.path.splitext(inputfile)[0] + "_headers.txt"
#name of the fasta file that doesn't have \n lines
stripped_fa = "stripped_" + inputfile2
#name of the final output fasta file
out_fasta = "filtered_" + inputfile2

#determine the clusters and subset the cd-hit file for only a few swine and human samples.
with open(inputfile) as infile, open (outputfile, "w") as outfile:
    lines = infile.readlines()
    cluster_lines = []
    for iter, line in enumerate(lines):
        if ">Cluster" in line:
            cluster_lines.append(iter)
    cluster_lines.append(len(lines))

    for i in range(len(cluster_lines)-1):
        swine_list = []
        human_list = []
        swine = 0
        human = 0
        for line in lines[cluster_lines[i]:cluster_lines[i+1]]:
            if ">Cluster" in line:
                pass
            elif ("/swine/") in line:
                swine_list.append(line)
                swine += 1
            elif ("/Swine/") in line:
                swine_list.append(line)
                swine += 1
            elif ("/sw/") in line:
                swine_list.append(line)
                swine += 1
            elif ("/SW/") in line:
                swine_list.append(line)
                swine += 1
            elif ("/sus_scrofa/") in line:
                swine_list.append(line)
                swine += 1
            elif ("/Sus_scrofa/") in line:
                swine_list.append(line)
                swine += 1
            elif ("/pig/") in line:
                swine_list.append(line)
                swine += 1
            elif ("/Pig/") in line:
                swine_list.append(line)
                swine += 1
            elif ("wild_boar") in line:
                swine_list.append(line)
                swine += 1
            else:
                human_list.append(line)
                human +=1
        outfile.write(">Cluster " + str(i) + "\n")
        if swine > 3 and human >3:
            for i in range(3):
                outfile.write(str(human_list[i]))
                outfile.write(str(swine_list[i]))
        elif 1 <= swine <=3 and 1 <= human <=3:
            for i in range(human):
                outfile.write(str(human_list[i]))
            for i in range(swine):
                outfile.write(str(swine_list[i]))
        elif 1 <= swine <=3 and human >3:
            for i in range(3):
                outfile.write(str(human_list[i]))
            for i in range(swine):
                outfile.write(str(swine_list[i]))
        elif swine > 3 and 1 <= human <=3:
            for i in range(human):
                outfile.write(str(human_list[i]))
            for i in range(3):
                outfile.write(str(swine_list[i]))
        elif swine >= 1 and human == 0:
            outfile.write(str(swine_list[0]))
        elif human >= 1 and swine == 0:
            outfile.write(str(human_list[0]))

#gather only the name of the sequences included in the subset. 
with open(outputfile) as out, open(outputfile2, "w") as out2:
    lines = out.readlines()
    for line in lines:
        if ">Cluster" in line:
            pass
        elif ">Cluster" not in line:
            if ">" in line:
                id1 = line.index(">")
                id2 = line.index("...")
                res = line[id1 + len("") + 1: id2]
                out2.write(">" + res + "\n")
        else:
            pass

#remove extra lines from the original fasta file
with open(inputfile2, "r") as in2, open(stripped_fa, "w") as out:
    lines = in2.readlines()
    for line in lines:
        if ">" not in line:
            out.write(line.strip())
        else:
            out.write("\n" + line)
            
#generate a new fasta file having only the subset from the cd-hit + the extra added sequences
with open(stripped_fa, "r") as in_fa, open(outputfile2, "r") as headers, open(out_fasta, "w") as out_fasta :
    lines = in_fa.readlines()
    lines_headers = headers.readlines()
    for index, line in enumerate(lines):
        for h_lines in lines_headers:
            if h_lines in line:
                out_fasta.write("".join(lines[max(0, index):index + 2])) 