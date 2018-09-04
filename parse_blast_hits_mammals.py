# Modules
from __future__ import division
import sys
# Input parameters
blast_file = sys.argv[1]
mam_pcid = int(sys.argv[2]) # Minimum percentage of identities for mammals genes
q_clen = int(sys.argv[3]) # Minimum percentage of matching cumulative length for the query
output_file = sys.argv[4]
# Functions
def FilterByPercentageOfIdentities(line):
    parse = line[:-1].split("\t")
    pcid = float(parse[2])
    # Fixed threshold for percentage of identities
    if pcid >= mam_pcid:
        test = True
    else:
        test = False
    return test
def StoreHitsPerQueryPerSubject(line, hits):
    parse = line[:-1].split("\t")
    query, subject = parse[0:2]
    # Parse subject accession IDs
    sallseqid = parse[15].split(";")
    prefix = [seqid.split("|")[1][0:3] if "|" in seqid else seqid[0:3] for seqid in sallseqid]
    # Accession prefixes for protein molecules
    if "NP_" in prefix or "XP_" in prefix or "AP_" in prefix or "YP_" in prefix or "WP_" in prefix:
        if not hits.has_key(query):
            hits[query] = {}
        if not hits[query].has_key(subject):
            hits[query][subject] = []
        hits[query][subject].append(parse)
    return query, hits
def FilterByCumulativeLength(query, hits, q_clen):
    query_matches = [[], []]
    for subject in hits[query]:
        q_len = int(hits[query][subject][0][12])
        q_cov = [0] * q_len
        for hit in hits[query][subject]:
            q_start, q_stop = map(int, hit[6:8])
            # Parse query coverage
            if 0 in q_cov[q_start-1:q_stop]:
                for i in range(q_start-1, q_stop):
                    q_cov[i] = 1
        # Query cumulative length
        q_tot_cov = len([i for i, x in enumerate(q_cov) if x == 1])
        q_match = (q_tot_cov/q_len) * 100
        # Fixed threshold for query cumulative length
        if q_match >= q_clen:
            sallseqid = hit[15].split(";")
            staxids =  hit[17].split(";")
            # Human genes
            if staxids == ["9606"]:
                for seqid in sallseqid:
                    query_matches[0].append(seqid.split("|")[1] if "|" in seqid else seqid)
            # Mouse genes
            elif staxids == ["10090"]:
                for seqid in sallseqid:
                    query_matches[1].append(seqid.split("|")[1] if "|" in seqid else seqid)
            # Human and mouse genes
            else:
                salltitles = hit[16].split("<>")
                for i in range(len(sallseqid)):
                    seqid = sallseqid[i]
                    if "Homo sapiens" in salltitles[i]:
                        query_matches[0].append(seqid.split("|")[1] if "|" in seqid else seqid)
                    elif "Mus musculus" in salltitles[i]:
                        query_matches[1].append(seqid.split("|")[1] if "|" in seqid else seqid)
    return query_matches
# Parse BLAST hits
hits = {}
file_out = open(output_file, "w")
file_out.write("# GeneID\tHuman\tMouse\n")
for line in open(blast_file, "rU"):
    if line[0] != "#":
        test = True
        # Filter hits by percentage of identities
        if test:
            test = FilterByPercentageOfIdentities(line)
        # Store hits per query sequence per subject sequence
        if test:
            query, hits = StoreHitsPerQueryPerSubject(line, hits)
    elif line[0] == "#" and hits != {}:
        # Filter hits per cumulative length of the query
        query_matches = FilterByCumulativeLength(query, hits, q_clen)
        # Output filtered hits
        file_out.write("\t".join([query, ";".join(query_matches[0]) if query_matches[0] != [] else "*", ";".join(query_matches[1]) if query_matches[1] != [] else "*"])+"\n")
        # Re-initialize variable
        hits = {}
file_out.close()
