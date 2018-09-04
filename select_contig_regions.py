# Modules
import sys
# Input parameters
fasta_file, assigned_file, unassigned_file, unmapped_file = sys.argv[1:5]
min_length = int(sys.argv[5]) # Minimum length of continuous uncovered base pairs
output_file = sys.argv[6]
# Extract contig sequences
all_contigs = {}
all_sequences = {}
for line in open(fasta_file, "rU"):
    if line[0] == ">":
        transcript = line[1:-1].split(" ")[0]
        gene = "_".join(transcript.split("_")[0:2])
        # Store gene and transcript IDs
        if not all_contigs.has_key(gene):
            all_contigs[gene] = []
        all_contigs[gene].append(transcript)
        all_sequences[transcript] = ""
    else:
        all_sequences[transcript] += line[:-1]
# Parse processed contigs
parsed = {}
# Assigned contigs
for line in open(assigned_file, "rU"):
    if line[0] != "#":
        transcript, hit, regions = line[:-1].split("\t")
        gene = "_".join(transcript.split("_")[0:2])
        if not parsed.has_key(gene):
            parsed[gene] = {}
        # Contig fully covered
        if regions == "0":
            parsed[gene][transcript] = ["0"]
        # Contig partly covered
        else:
            match = []
            for region in regions.split(";"):
                start, stop = map(int, region.split("-"))
                # At least "min_length" continuous base pairs must be uncovered
                if (stop-start+1) >= min_length:
                    match.append(region)
            if len(match) == 0:
                match = ["0"]
            parsed[gene][transcript] = match
# Unassigned contigs
for line in open(unassigned_file, "rU"):
    if line[0] != "#":
        transcript, hit, regions = line[:-1].split("\t")
        gene = "_".join(transcript.split("_")[0:2])
        if not parsed.has_key(gene):
            parsed[gene] = {}
        parsed[gene][transcript] = [regions]
# Unmapped contigs
for line in open(unmapped_file, "rU"):
    if line[0] != "#":
        transcript, hit, regions = line[:-1].split("\t")
        gene = "_".join(transcript.split("_")[0:2])
        if not parsed.has_key(gene):
            parsed[gene] = {}
        parsed[gene][transcript] = [regions]
# Output sequences of contig regions
file_out = open(output_file, "w")
for gene in sorted(parsed):
    for transcript in sorted(all_contigs[gene]):
        headers = []
        seqs = []
        # Unparsed, unassigned and unmapped contigs 
        if not parsed[gene].has_key(transcript) or parsed[gene][transcript] == ["*"]:
            seqs = [all_sequences[transcript]]
            headers = [">"+gene+"|"+transcript+" "+"len="+str(len(seqs[0]))]
        # Assigned contigs
        else:
            inc = 0
            for region in parsed[gene][transcript]:
                # Uncovered regions
                if region != "0":
                    inc += 1
                    start, stop = map(int, region.split("-"))
                    seq = all_sequences[transcript][start-1:stop]
                    seqs.append(seq)
                    header = ">"+gene+"|"+transcript+"|"+"r"+str(inc)+"_"+region+" "+"len="+str(len(seq))
                    headers.append(header)
        # Output headers and sequences
        for i in range(len(headers)):
            # Uncovered regions of at least "min_length" base pairs
            if len(seqs[i]) >= min_length:
                output = headers[i]+"\n"
                for b in range((len(seqs[i])/60)+1):
                    text = seqs[i][b*60:(b+1)*60]
                    output += text
                    if text != "":
                        output += "\n"
                file_out.write(output)
file_out.close()
