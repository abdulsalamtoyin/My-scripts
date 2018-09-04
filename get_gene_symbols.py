### Commands arising from IITA bioinformatics Resources Note, by Abdulsalam Toyin.
# Modules
import sys
# Input parameters
taxids_file, refseq_file, db_type = sys.argv[1:4]
hits_file = sys.argv[4]
output_file = sys.argv[5]
# Get taxids
taxids = {}
for line in open(taxids_file, "rU"):
    taxids[line[:-1]] = ""
# Get gene symbols
symbols = {}
for line in open(refseq_file, "rU"):
    d = line[:-1].split("\t")
    taxid = d[0]
    if db_type == "nucl":
        accession = d[3]
    elif db_type == "prot":
        accession = d[5]
    symbol = d[15]
    if taxids.has_key(taxid) and accession != "-" and symbol != "-":
        if not symbols.has_key(accession):
            symbols[accession] = [symbol]
        elif not symbol in symbols[accession]:
            symbols[accession].append(symbol)
# Assign symbols to hits
file_out = open(output_file, "w")
for line in open(hits_file, "rU"):
    if line[0] != "#":
        gene = line[:-1].split("\t")[0]
        hits = line[:-1].split("\t")[1:]
        l_symbols = [[] for i in range(len(hits))]
        # Parse hits
        for i in range(len(hits)):
            if hits[i] != "*":
                for hit in hits[i].split(";"):
                    # Parse symbols
                    if symbols.has_key(hit):
                        for symbol in symbols[hit]:
                            if not symbol in l_symbols[i]:
                                l_symbols[i].append(symbol)
        # Output gene symbol assignment
        file_out.write(gene+"\t"+"\t".join([";".join(l) if l != [] else "*" for l in l_symbols])+"\n")
    else:
        file_out.write(line)
file_out.close()
