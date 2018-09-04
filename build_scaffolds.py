### Commands arising from IITA bioinformatics Resources Note, by Abdulsalam Toyin.
# Modules
from itertools import groupby
from operator import itemgetter
import ast
import sys
# Input parameters
index_file, blast_file, output_file = sys.argv[1:]
# Load index
g_index = {}
for line in open(index_file, "rU"):
    if line[0] != "#":
        d = line[:-1].split("\t")
        gene = d[0]
        nb_regions = int(d[1])
        isoform_ids = d[2].split(";")
        region_ids = d[3].split(";")
        regions = []
        for region in d[4:]:
            regions.append(ast.literal_eval(region))
        g_index[gene] = [nb_regions, isoform_ids, region_ids, regions]
# Load hits from BLAST file
to_parse = {}
scaffolds = {}
hits = {}
q_length = {}
for line in open(blast_file, "rU"):
    if line[0] != "#":
        d = line[:-1].split("\t")
        q_gene = d[0].split("|")[0]
        q_isoform = "|".join(d[0].split("|")[1:])
        s_gene = d[1].split("|")[0]
        s_isoform = "|".join(d[1].split("|")[1:])
        q_len, s_len = map(int, d[12:14])
        # Store the length of the query
        q_length[q_isoform] = q_len
        # Gene isoform according to Trinity
        if q_gene == s_gene:
            # Isoform that has not been parsed yet
            if not q_isoform in g_index[q_gene][2]:
                # Extract isoform that is entirely matched with a single hit
                q_start, q_stop = map(int, d[6:8])
                if q_start == 1 and q_stop == q_len:
                    if not scaffolds.has_key(q_gene):
                        scaffolds[q_gene] = {}
                    if not scaffolds[q_gene].has_key(s_isoform):
                        scaffolds[q_gene][s_isoform] = {}
                    if not scaffolds[q_gene][s_isoform].has_key(q_isoform):
                        scaffolds[q_gene][s_isoform][q_isoform] = [[[s_isoform, 1, s_len], [q_isoform, 1, q_len]]]
                # Isoform that needs to be parsed
                elif not scaffolds.has_key(q_gene) or not scaffolds[q_gene].has_key(s_isoform) or not scaffolds[q_gene][s_isoform].has_key(q_isoform):
                    q_transcript_id = q_isoform.split("|")[0]
                    # Store hit(s) for the selected isoform
                    if not to_parse.has_key(q_gene):
                        to_parse[q_gene] = q_transcript_id
                        hits[q_gene] = [d]
                    elif to_parse[q_gene] == q_transcript_id:
                        hits[q_gene].append(d)
# Parse hits
coverage = {}
for gene in to_parse:
    for hit in hits[gene]:
        q_isoform = "|".join(hit[0].split("|")[1:])
        s_isoform = "|".join(hit[1].split("|")[1:])
        # Store gene
        if not coverage.has_key(gene):
            coverage[gene] = {}
        # Store subject isoform per region
        if not coverage[gene].has_key(s_isoform):
            coverage[gene][s_isoform] = {}
        # Initialize subject isoform coverage per query isoform
        if not coverage[gene][s_isoform].has_key(q_isoform):
            s_len = int(hit[13])
            coverage[gene][s_isoform][q_isoform] = [0] * s_len
        # Coverage of the query per subject sequence
        q_pos = int(hit[6])
        s_start, s_stop = map(int, hit[8:10])
        for s_pos in range(s_start-1, s_stop):
            if coverage[gene][s_isoform][q_isoform][s_pos] == 0:
                coverage[gene][s_isoform][q_isoform][s_pos] = q_pos
            q_pos += 1
# Parse coverage
for gene in coverage:
    if not scaffolds.has_key(gene):
        scaffolds[gene] = {}
    for s_isoform in coverage[gene]:
        if not scaffolds[gene].has_key(s_isoform):
            scaffolds[gene][s_isoform] = {}
        for q_isoform in coverage[gene][s_isoform]:
            consensus = []
            # Positions of the query sequence
            q_matched = filter(None, coverage[gene][s_isoform][q_isoform])
            set_query = set(q_matched)
            q_unmatched = [pos for pos in range(1, q_length[q_isoform]+1) if pos not in set_query]
            # The query sequence is entirely matched
            if q_unmatched == []:
                consensus.append([[s_isoform, 1, len(coverage[gene][s_isoform][q_isoform])], [q_isoform, 1, q_length[q_isoform]]])
            # The query sequence is partly unmatched
            else:
                # Identify the unmatched regions of the query sequence 
                q_unmatched_regions = []
                for k, g in groupby(enumerate(q_unmatched), lambda (i, x): i-x):
                    pos = map(itemgetter(1), g) 
                    q_unmatched_regions.append(pos[0])
                    q_unmatched_regions.append(pos[-1])
                # The beginning of the query sequence does not match
                q_start, q_stop = 1, min(q_matched)
                if q_stop != 1:
                    consensus.append([[], [q_isoform, q_start, q_stop-1]])
                # Parse the coverage of the subject sequence
                s_length = len(coverage[gene][s_isoform][q_isoform])
                pos = 0
                while pos < s_length:
                    q_start = coverage[gene][s_isoform][q_isoform][pos]
                    s_start = pos
                    # The subject sequence is not covered
                    if q_start == 0:
                        while pos < s_length and coverage[gene][s_isoform][q_isoform][pos] == 0:
                            pos += 1
                        s_stop = pos
                        consensus.append([[s_isoform, s_start+1, s_stop], []])
                    else:
                        # Unmatched region of the query sequence
                        if q_start != q_stop:
                            if (q_stop+1) in q_unmatched_regions and (q_start-1) in q_unmatched_regions:
                                consensus.append([[], [q_isoform, q_stop+1, q_start-1]])
                        # The query sequence matches the subject sequence
                        while pos < s_length and coverage[gene][s_isoform][q_isoform][pos] != 0:
                            # Unmatched gap within the query sequence
                            gap_start = coverage[gene][s_isoform][q_isoform][pos]+1
                            if (pos+1) < s_length:
                                gap_stop = coverage[gene][s_isoform][q_isoform][pos+1]-1 
                                if gap_start in q_unmatched_regions and gap_stop in q_unmatched_regions:
                                    s_stop = pos+1
                                    q_stop = coverage[gene][s_isoform][q_isoform][pos]
                                    consensus.append([[s_isoform, s_start+1, s_stop], [q_isoform, q_start, q_stop]])
                                    consensus.append([[], [q_isoform, gap_start, gap_stop]])
                                    q_start = coverage[gene][s_isoform][q_isoform][pos+1]
                                    s_start = pos+1
                            pos += 1
                        s_stop = pos
                        q_stop = coverage[gene][s_isoform][q_isoform][pos-1]
                        consensus.append([[s_isoform, s_start+1, s_stop], [q_isoform, q_start, q_stop]])
                # The end of the query sequence does not match
                q_stop = max(q_matched)
                if q_stop != q_length[q_isoform]:
                    consensus.append([[], [q_isoform, q_stop+1, q_length[q_isoform]]])
            # Store consensus coordinates
            if not scaffolds[gene][s_isoform].has_key(q_isoform):
                scaffolds[gene][s_isoform][q_isoform] = consensus
# Export scaffolds
file_out = open(output_file, "w")
file_out.write("# TrinityID\tSubjectID\tQueryID\tScaffold\n")
for gene in sorted(scaffolds):
    for s_isoform in sorted(scaffolds[gene]):
        for q_isoform in sorted(scaffolds[gene][s_isoform]):
            file_out.write("\t".join([gene, s_isoform, q_isoform, str(scaffolds[gene][s_isoform][q_isoform])])+"\n")
file_out.close()
