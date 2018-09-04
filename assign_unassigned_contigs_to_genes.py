### Commands arising from IITA bioinformatics Resources Note, by Abdulsalam Toyin.
# Modules
import sys
# Input parameters
bed_file, output_file = sys.argv[1:]
# Parse BED file
hits = {}
for line in open(bed_file, "rU"):
    d = line[:-1].split("\t")
    trinity_id = d[3]
    matched_id = d[9]
    if not hits.has_key(trinity_id):
        hits[trinity_id] = []
    if not matched_id in hits[trinity_id]:
        hits[trinity_id].append(matched_id)
# Output assignment
file_out = open(output_file, "w")
file_out.write("# TrinityID\tGeneID\tUnmatchedRegions\n")
for trinity_id in sorted(hits):
    file_out.write("\t".join([trinity_id, ";".join(hits[trinity_id]), "*"])+"\n")
file_out.close()
