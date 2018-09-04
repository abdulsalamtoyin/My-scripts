### Commands arising from IITA bioinformatics Resources Note, by Abdulsalam Toyin.
###2017

# Modules
import sys
# Input parameters
gtf_file, output_file = sys.argv[1:]
# Parse GTF file
genes = {}
for line in open(gtf_file, "rU"):
    d = line[:-2].split("\t")
    # Extract exon coordinates
    chr = d[0]
    start, stop = map(int, d[3:5])
    strand = d[6]
    ids = d[8].split("; ")
    for id in ids:
        if id[0:7] == "gene_id":
            gid = id[9:-1]
    # Store exon coordinates
    if not genes.has_key(gid):
        genes[gid] = {}
    # Gene IDs spanning on multiple chromosomes
    if not genes[gid].has_key(chr):
        genes[gid][chr] = [strand, [start, stop]]
    else:
        genes[gid][chr][1] += [start, stop]
# Output gene coordinates
file_out = open(output_file, "w")
for gid in sorted(genes):
    for chr in sorted(genes[gid]):
        strand = genes[gid][chr][0]
        start = str(min(genes[gid][chr][1]))
        stop = str(max(genes[gid][chr][1]))
        file_out.write("\t".join([chr, start, stop, gid, ".", strand])+"\n")
file_out.close()
