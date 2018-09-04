### Commands arising from IITA bioinformatics Resources Note, by Abdulsalam Toyin.
# Modules
import sys
# Input parameters
fasta_file, output_file = sys.argv[1:]
# Parse gene and transcript IDs
file_out = open(output_file, "w")
for line in open(fasta_file, "rU"):
    if line[0] == ">":
        transcript = line[1:-1].split(" ")[0]
        gene = transcript.split("|")[0]
        file_out.write("\t".join([gene, transcript])+"\n")
file_out.close()
