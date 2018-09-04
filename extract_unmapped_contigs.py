### Commands arising from IITA bioinformatics Resources Note, by Abdulsalam Toyin.
# Modules
import sys
# Input parameters
fasta_file, assigned_file, unassigned_file, output_file = sys.argv[1:]
# Extract assigned genes
assigned = {}
for line in open(assigned_file, "rU"):
    if line[0] != "#":
        gene = "_".join(line[:-1].split("\t")[0].split("_")[0:2])
        assigned[gene] = ""
# Extract assigned genes
unassigned = {}
for line in open(unassigned_file, "rU"):
    if line[0] != "#":
        gene = "_".join(line[:-1].split("\t")[0].split("_")[0:2])
        unassigned[gene] = ""
# Extract unmapped genes
file_out = open(output_file, "w")
file_out.write("# TrinityID\tGeneID\tUnmatchedRegions\n")
for line in open(fasta_file, "rU"):
    if line[0] == ">":
        transcript = line[1:-1].split(" ")[0]
        gene = "_".join(transcript.split("_")[0:2])
        # Genes with all isoforms unmapped
        if not assigned.has_key(gene) and not unassigned.has_key(gene):
            file_out.write("\t".join([transcript, "*", "*"])+"\n")
file_out.close()
