### Commands arising from IITA bioinformatics Resources Note, by Abdulsalam Toyin.
# Modules
from itertools import groupby
from operator import itemgetter
import sys
# Input parameters
blast_file = sys.argv[1]
t_aln_length = int(sys.argv[2]) # Minimum number of alignment base pairs
output_file = sys.argv[3]
# Functions
def FilterByAlignmentLength(line, t_aln_length):
    parse = line[:-1].split("\t")
    aln_length = int(parse[3])
    # Fixed threshold for alignment length
    if aln_length >= t_aln_length:
        test = True
    else:
        test = False
    return test
def ParseHitsPerQuery(line, t_aln_length, hits, query_cumul_length):
    parse = line[:-1].split("\t")
    query, subject = parse[0:2]
    start, stop = map(int, parse[6:8])
    test = False
    if not query_cumul_length.has_key(query):
        # Set query coverage
        query_cumul_length[query] = [1] * int(parse[12])
        test = True
        # Store subject ID
        hits[query] = [subject]
    else:
        # Unique subject ID matching
        if hits[query] == [subject]:
            if 1 in query_cumul_length[query][start-1:stop]:
                test = True
        # Multiple subject IDs matching
        else:
            # At least "t_aln_length" bp must be uncovered from previous hits 
            if 1 in query_cumul_length[query][start-1:stop] and sum(query_cumul_length[query][start-1:stop]) >= t_aln_length:
                test = True
                # Store subject ID
                if not subject in hits[query]:
                    hits[query].append(subject)
    # Update query coverage
    if test:
        for i in range(start-1, stop):
            query_cumul_length[query][i] = 0
    return hits, query_cumul_length
# Parse BLAST hits
hits = {}
query_cumul_length = {}
data = open(blast_file, "rU").readlines()
file_out = open(output_file, "w")
file_out.write("# TrinityID\tGeneID\tUnmatchedRegions\n")
for line in range(len(data)):
    if data[line][0] != "#":
        test = True
        # Filter hits by alignment length
        if test:
            test = FilterByAlignmentLength(data[line], int(t_aln_length))
        # Parse hits per query sequence
        if test:
            hits, query_cumul_length = ParseHitsPerQuery(data[line], t_aln_length, hits, query_cumul_length)
        # Output filtered hits
        query = data[line][:-1].split("\t")[0]
        if query_cumul_length.has_key(query) and (line == (len(data)-1) or data[line+1][0] == "#" or hits.has_key(data[line+1][:-1].split("\t")[0]) == False):
            file_out.write(query+"\t"+";".join(hits[query])+"\t")
            # Parse query coverage to retrieve unmatched positions
            unmatched = [i for i, x in enumerate(query_cumul_length[query]) if x == 1]
            # Query sequence is fully covered
            if unmatched == []:
                file_out.write("0"+"\n")
            # Query sequence has unmatched regions
            else:
                regions = []
                # Concatenate consecutive regions
                for k, g in groupby(enumerate(unmatched), lambda (i, x): i-x):
                    pos = map(itemgetter(1), g) 
                    regions.append(str(pos[0]+1)+"-"+str(pos[-1]+1))
                file_out.write(";".join(regions)+"\n")
            # Re-initialize variables
            hits = {}
            query_cumul_length = {}
file_out.close()
