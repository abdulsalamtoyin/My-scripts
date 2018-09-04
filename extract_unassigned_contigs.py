### Commands arising from IITA bioinformatics Resources Note, by Abdulsalam Toyin.
# Modules
import sys
# Input parameters
fasta_file, assigned_file, output_file = sys.argv[1:]
# Extract assigned genes
assigned = {}
for line in open(assigned_file, "rU"):
    if line[0] != "#":
        gene = "_".join(line[:-1].split("\t")[0].split("_")[0:2])
        assigned[gene] = ""
# Extract unassigned genes
test = False
file_out = open(output_file, "w")
for line in open(fasta_file, "rU"):
    if line[0] == ">":
        gene = "_".join(line[1:-1].split(" ")[0].split("_")[0:2])
        # Genes with all isoforms unassigned
        if not assigned.has_key(gene):
            file_out.write(line)
            test = True
        else:
            test = False
    elif test:
        file_out.write(line)
file_out.close()
