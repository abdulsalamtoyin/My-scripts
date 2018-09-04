### Commands arising from IITA bioinformatics Resources Note, by Abdulsalam Toyin.
# Modules
import sys
# Input parameters
gtf_file, blast_file, output_file = sys.argv[1:]
# Extract gene and transcript IDs
transcript_ids = {} 
for line in open(gtf_file, "rU"):
    ids = line[:-2].split("\t")[8].split("; ")
    for id in ids:
        if id[0:7] == "gene_id":
            gid = id[9:-1]
        elif id[0:13] == "transcript_id":
            tid = id[15:-1]
	    #print(tid)
    transcript_ids[tid] = gid
    #print(transcript_ids[tid] = gid)
# Parse BLAST file
file_out = open(output_file, "w") 
for line in open(blast_file, "rU"):
    if line[0] == "#":
        file_out.write(line)
    else:
        d = line.split("\t")
	#print(d)
        #d[1] = transcript_ids[d[1]]
        file_out.write("\t".join(d)) 
file_out.close()
