### Commands arising from IITA bioinformatics Resources Note, by Abdulsalam Toyin.
# Modules
import os
import sys
# Input parameters
query_name, blast_name, filename = sys.argv[1:]
# Function
def CheckProcessing(file_in):
    nb_genes = processed = remaining = 0
    for line in open(file_in, "rU"):
        if line[0] != "#":
            nb_genes += 1
            d = line[:-1].split("\t")
            nb_regions = int(d[1])
            regions_ids = d[3].split(";")
            # Procesed genes
            if len(regions_ids) == nb_regions:
                processed += 1
            # Remaining genes
            else:
                remaining += 1
    print "= Processed: %i/%i | Remaining: %i/%i\n" % (processed, nb_genes, remaining, nb_genes)
    return remaining
# Incremental gene concatenation
processing = True
inc = 0
while processing == True:
    inc += 1
    prefix = str(inc)+"_"
    print "* Processing: #"+str(inc)
    # Create a BLAST database: makeblastdb
    print "> Create a BLAST database: makeblastdb"
    cmd = ["makeblastdb", "-in", prefix+filename+".fa", "-dbtype", "nucl", "-logfile", "makeblastdb.log"]
    os.system(" ".join(cmd))
    # Compare sequences: blastn
    print "> Compare contig sequences: blastn"
    options = " ".join(["-num_threads 12", "-perc_identity", "90", "-strand", "plus", "-dust", "no", "-soft_masking", "no", "-ungapped", "-outfmt", "\"7 std qlen slen sstrand\""])
    cmd = ["blastn", "-query", query_name, "-db", prefix+filename+".fa", options , "-out", prefix+blast_name]
    os.system(" ".join(cmd))
    # Build contig scaffolds
    print "> Build contig scaffolds"
    cmd = ["python", "script/build_scaffolds.py", prefix+filename+".index", prefix+blast_name, prefix+filename+".scaffolds"]
    os.system(" ".join(cmd))
    # Build contig sequences
    print "> Build contig sequences"
    inname = prefix+filename
    outname = str(inc+1)+"_"+filename
    cmd = ["python", "script/build_consensus_sequences.py", query_name, inname+".scaffolds", inname+".fa", inname+".index", outname+".fa", outname+".index"]
    os.system(" ".join(cmd))
    # Check processing
    remaining = CheckProcessing(outname+".index")
    if remaining == 0:
        processing = False
   # break
        print "Final output files:\n%s\n%s" % (outname+".index", outname+".fa")
    #break
    #print("finish")
