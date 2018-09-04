# Modules
from __future__ import division
import sys
# Input parameters
blast_file = sys.argv[1]
ggal_pcid = int(sys.argv[2]) # Minimum percentage of identities for chicken genes (staxid == 9031)
other_pcid = int(sys.argv[3]) # Minimum percentage of identities for other bird genes (staxid != 9031)
q_clen = int(sys.argv[4]) # Minimum percentage of matching cumulative length for the query
s_clen = int(sys.argv[5]) # Minimum percentage of matching cumulative length for the subject
output_file = sys.argv[6]
# Functions
def FilterByPercentageOfIdentities(line):
    parse = line[:-1].split("\t")
    pcid = float(parse[2])
    staxids = parse[17].split(";")
    # Fixed threshold for percentage of identities
    if "9031" in staxids and pcid >= ggal_pcid:
        test = True
    elif staxids != ["9031"] and pcid >= other_pcid:
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
    # Accession prefixes for RNA molecules
    if "NM_" in prefix or "NR_" in prefix or "XM_" in prefix or "XR_" in prefix:
        if not hits.has_key(query):
            hits[query] = {}
        if not hits[query].has_key(subject):
            hits[query][subject] = []
        hits[query][subject].append(parse)
    return query, hits
def FilterByCumulativeLength(query, hits, q_clen, s_clen):
    query_matches = [[], []]
    for subject in hits[query]:
        q_len, s_len = map(int, hits[query][subject][0][12:14])
        q_cov = [0] * q_len
        s_cov = [0] * s_len
        for hit in hits[query][subject]:
            q_start, q_stop, s_start, s_stop = map(int, hit[6:10])
            # Parse query coverage
            if 0 in q_cov[q_start-1:q_stop]:
                for i in range(q_start-1, q_stop):
                    q_cov[i] = 1
            # Parse subject coverage
            if 0 in s_cov[s_start-1:s_stop]:
                for i in range(s_start-1, s_stop):
                    s_cov[i] = 1
        # Query cumulative length
        q_tot_cov = len([i for i, x in enumerate(q_cov) if x == 1])
        q_match = (q_tot_cov/q_len) * 100
        # Subject cumulative length
        s_tot_cov = len([i for i, x in enumerate(s_cov) if x == 1])
        s_match = (s_tot_cov/s_len) * 100
        # Fixed threshold for query/subject cumulative length
        if q_match >= q_clen and s_match >= s_clen:
            sallseqid = hit[15].split(";")
            staxids =  hit[17].split(";")
            # Chicken genes
            if staxids == ["9031"]:
                for seqid in sallseqid:
                    query_matches[0].append(seqid.split("|")[1] if "|" in seqid else seqid)
            # Genes from chicken and other birds
            elif "9031" in staxids:
                salltitles = hit[16].split("<>")
                for i in range(len(sallseqid)):
                    seqid = sallseqid[i]
                    if "Gallus gallus" in salltitles[i]:
                        query_matches[0].append(seqid.split("|")[1] if "|" in seqid else seqid)
                    else:
                        query_matches[1].append(seqid.split("|")[1] if "|" in seqid else seqid)
            # Genes from other birds
            else:
                for seqid in sallseqid:
                    query_matches[1].append(seqid.split("|")[1] if "|" in seqid else seqid)
    return query_matches
# Parse BLAST hits
hits = {}
file_out = open(output_file, "w")
file_out.write("# GeneID\tChicken\tBirds\n")
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
        # Filter hits per cumulative length of the query and the subject
        query_matches = FilterByCumulativeLength(query, hits, q_clen, s_clen)
        # Output filtered hits
        file_out.write("\t".join([query, ";".join(query_matches[0]) if query_matches[0] != [] else "*", ";".join(query_matches[1]) if query_matches[1] != [] else "*"])+"\n")
        # Re-initialize variable
        hits = {}
file_out.close()
