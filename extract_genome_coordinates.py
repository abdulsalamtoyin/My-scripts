### Commands arising from IITA bioinformatics Resources Note, by Abdulsalam Toyin.
# Modules
import sys
# Input parameters
input_file, output_file = sys.argv[1:]
# Parse BLAST file
file_out = open(output_file, "w")
for line in open(input_file, "rU"):
    if line[0] != "#":
        d = line[:-1].split("\t")
        gid, chr = d[0:2]
        if d[14] == "plus":
            strand = "+"
            start, stop = d[8:10]
        elif d[14] == "minus":
            strand = "-"
            stop, start = d[8:10]
        file_out.write("\t".join([chr, start, stop, gid, ".", strand])+"\n")
file_out.close()
