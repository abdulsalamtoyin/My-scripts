### Commands arising from IITA bioinformatics Resources Note, by Abdulsalam Toyin.
# Modules
import sys
# Input parameters
map_file, db_file, output_file = sys.argv[1:]
# Get seqids
seqids = {}
for line in open(map_file, "rU"):
    seqid, taxid = line[:-1].split("\t")
    seqids[seqid] = taxid
# Extract database sequences
file_out = open(output_file, "w")
for line in open(db_file, "rU"):
    if line[0] == ">":
        test = False
        headers = []
        for header in line[1:-1].split("\x01"):
            seqid = header.split(" ")[0]
            if seqids.has_key(seqid):
                test = True
                headers.append(header)
        # Headers
        if test == True:
            file_out.write(">"+"\x01".join(headers)+"\n")
    # Sequences
    elif test == True:
        file_out.write(line)
file_out.close()
