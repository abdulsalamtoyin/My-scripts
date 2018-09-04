### Commands arising from IITA bioinformatics Resources Note, by Abdulsalam Toyin.
# Modules
import sys
# Input parameters
taxids_file, map_file, output_file = sys.argv[1:]
# Get taxids
taxids = {}
for line in open(taxids_file, "rU"):
    taxids[line[:-1]] = ""
# Extract accession IDs
file_out = open(output_file, "w")
for line in open(map_file, "rU"):
    accession, taxid = line[:-1].split("\t")[1:3] 
    if taxids.has_key(taxid):
        file_out.write("\t".join([accession, taxid])+"\n")
file_out.close()
