### Commands arising from IITA bioinformatics Resources Note, by Abdulsalam Toyin.
# Modules
import sys
# Input parameters
assignment_file, fasta_file, output_file = sys.argv[1:]
# Parse assignment
assignment = {}
for line in open(assignment_file, "rU"):
    if line[0] != "#":
        gene = line[:-1].split("\t")[0]
        hits = line[:-1].split("\t")[1:]
        for hit in hits:
            if hit != "*":
                assignment[gene] = ""
                break
# Extract unassigned genes
file_out = open(output_file, "w")
for line in open(fasta_file, "rU"):
    if line[0] == ">":
        header = ""
        gene = line[1:-1].split(" ")[0]
        if not assignment.has_key(gene):
            header = line
            file_out.write(header)
    elif header != "":
        file_out.write(line)
file_out.close()
