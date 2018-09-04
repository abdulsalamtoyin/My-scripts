### Commands arising from IITA bioinformatics Resources Note, by Abdulsalam Toyin.
# Modules
import sys
# Input parameters
fasta_file, output_file, index_file = sys.argv[1:]
# Parse FASTA file
isoforms = {}
sequences = {}
lengths = {}
for line in open(fasta_file, "rU"):
    # Header
    if line[0] == ">":
        header = line
        gene, isoform = line[1:-1].split(" ")[0].split("|")[0:2]
        length = int(line[1:-1].split(" ")[1].split("=")[1])
        # Store gene
        if not isoforms.has_key(gene):
            isoforms[gene] = []
        # Store isoform per gene
        if not isoform in isoforms[gene]:
            isoforms[gene].append(isoform)
            lengths[isoform] = []
        # Store length per isoform
        lengths[isoform].append(length)
        # Store header
        if not sequences.has_key(isoform):
            sequences[isoform] = {}
        sequences[isoform][header] = ""
    # Sequence
    else:
        sequences[isoform][header] += line
# Identify longest isoform per gene
longest_isoforms = {}
for gene in isoforms:
    longest = ["", 0]
    for isoform in isoforms[gene]:
        if sum(lengths[isoform]) > longest[1]:
            longest = [isoform, sum(lengths[isoform])]
    # Store longest isoform
    longest_isoforms[longest[0]] = longest[1]
# Output longest isoform per gene
file_out = open(output_file, "w")
for isoform in sorted(longest_isoforms):
    # Isoform cut into several regions
    for header in sorted(sequences[isoform]):
        sequence = sequences[isoform][header]
        file_out.write(header+sequence)
file_out.close()
# Create an index file
index_out = open(index_file, "w")
index_out.write("# TrinityID\tNbRegions\tIsoformsID\tRegionsID\tRegions\n")
for gene in sorted(isoforms):
    # Number of regions per gene
    nb_regions = len([i for isoform in isoforms[gene] for i in lengths[isoform]])
    # Longest isoform per gene
    for isoform in isoforms[gene]:
        if longest_isoforms.has_key(isoform):
            longest_isoform = isoform
    # Regions
    region_ids = []
    regions = []
    for header in sorted(sequences[longest_isoform]):
        d = header[1:-1].split(" ")
        region = "|".join(d[0].split("|")[1:])
        region_ids.append(region)
        length = int(d[1].split("=")[1])
        regions.append([region, 1, length])
    # Output index file
    index_out.write("\t".join([gene, str(nb_regions), longest_isoform, ";".join(region_ids), "\t".join(str(region) for region in regions)])+"\n")
index_out.close()
