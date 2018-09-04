### Commands arising from IITA bioinformatics Resources Note, by Abdulsalam Toyin.
# Modules
import sys
# Input parameters
fasta_file, single_file, multi_file = sys.argv[1:]
# Parse FASTA file
isoforms = {}
sequences = {}
for line in open(fasta_file, "rU"):
    # Header
    if line[0] == ">":
        header = line
        gene, isoform = line[1:-1].split(" ")[0].split("|")[0:2]
        # Store gene
        if not isoforms.has_key(gene):
            isoforms[gene] = []
        # Store isoform per gene
        if not isoform in isoforms[gene]:
            isoforms[gene].append(isoform)
        # Store header
        if not sequences.has_key(isoform):
            sequences[isoform] = {}
        sequences[isoform][header] = ""
    # Sequence
    else:
        sequences[isoform][header] += line
# Split genes according to the number of isoforms
single_out = open(single_file, "w")
multi_out = open(multi_file, "w")
for gene in sorted(isoforms):
    # Genes with a single isoform
    if len(isoforms[gene]) == 1:
        isoform = isoforms[gene][0]
        # Isoform cut into several regions
        for header in sorted(sequences[isoform]):
            sequence = sequences[isoform][header]
            single_out.write(header+sequence)
    # Genes with multiple isoforms
    else:
        for isoform in isoforms[gene]:
            # Isoform cut into several regions
            for header in sorted(sequences[isoform]):
                sequence = sequences[isoform][header]
                multi_out.write(header+sequence)
single_out.close()
multi_out.close()
