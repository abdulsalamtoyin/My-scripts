# Modules
from __future__ import division
import sys
# Input parameters
blast_file = sys.argv[1]
t_aln_length = int(sys.argv[2]) # Minimum number of overlapping base pairs
p_cumul_length = int(sys.argv[3]) # Minimum percentage of cumulative alignment length
output_file = sys.argv[4]
# Function
def FilterByCumulativeLengthOfTheQuery(line, query_cumul_length, t_aln_length):
    parse = line[:-1].split("\t")
    query = parse[0]
    start, stop = map(int, parse[6:8])
    if not query_cumul_length.has_key(query):
        # Set query coverage
        query_cumul_length[query] = [1] * int(parse[12])
    # At least "t_aln_length" bp must be uncovered from previous hits
    if 1 in query_cumul_length[query][start-1:stop] and sum(query_cumul_length[query][start-1:stop]) >= t_aln_length:
        # Update query coverage
        for i in range(start-1, stop):
            query_cumul_length[query][i] = 0
        test = True
    else:
        test = False
    return test, query_cumul_length
# Parse BLAST hits
hits = []
query_cumul_length = {}
data = open(blast_file, "rU").readlines()
file_out = open(output_file, "w")
for line in range(len(data)):
    if data[line][0] != "#":
        # Filter by cumulative length of the query
        test, query_cumul_length = FilterByCumulativeLengthOfTheQuery(data[line], query_cumul_length, t_aln_length)
        if test:
            hits.append(data[line])
        # Output filtered hits
        query = data[line][:-1].split("\t")[0]
        if query_cumul_length.has_key(query) and (line == (len(data)-1) or data[line+1][0] == "#" or query_cumul_length.has_key(data[line+1][:-1].split("\t")[0]) == False):
            # Parse query coverage to retrieve matched positions
            matched = [i for i, x in enumerate(query_cumul_length[query]) if x == 0]
            # Output hits passing fixed threshold of cumulative length
            if (len(matched)/len(query_cumul_length[query])*100) >= p_cumul_length:
                for hit in hits:
                    file_out.write(hit)
            # Re-initialize variables
            hits = []
            query_cumul_length = {}
file_out.close()
