# Modules
import ast
import sys
# Input parameters
fasta_file, scaffold_file, db_file, index_file, output_fasta, output_index = sys.argv[1:]
# Load contig sequences
sequences = {}
for line in open(fasta_file, "rU"):
    if line[0] == ">":
        header = line[1:-1].split(" ")[0].split("|")
        gene, isoform = header[0:2]
        if not sequences.has_key(gene):
            sequences[gene] = {}
        if not sequences[gene].has_key(isoform):
            sequences[gene][isoform] = {}
        if len(header) == 2:
            region = isoform
        else:
            region = "|".join(header[1:])
        sequences[gene][isoform][region] = ""
    else:
        sequences[gene][isoform][region] += line[:-1]
# Load scaffolds
scaffolds = {}
for line in open(scaffold_file, "rU"):
    if line[0] != "#":
        d = line[:-1].split("\t")
        gene, s_isoform, q_isoform = d[0:3]
        scaffold = ast.literal_eval(d[3])
        if not scaffolds.has_key(gene):
            scaffolds[gene] = {}
        if not scaffolds[gene].has_key(s_isoform):
            scaffolds[gene][s_isoform] = {}
        if not scaffolds[gene][s_isoform].has_key(q_isoform):
            scaffolds[gene][s_isoform][q_isoform] = scaffold
# Load reference contig sequences
db = {}
for line in open(db_file, "rU"):
    if line[0] == ">":
        header = line[1:-1].split(" ")[0].split("|")
        gene = header[0]
        if not db.has_key(gene):
            db[gene] = {}
        scaffold = "|".join(header[1:])
        db[gene][scaffold] = ""
    else:
        db[gene][scaffold] += line[:-1]
# Build consensus sequences
fasta_out = open(output_fasta, "w")
index_out = open(output_index, "w")
for line in open(index_file, "rU"):
    if line[0] == "#":
        index_out.write(line)
    else:
        # Parse index
        d = line[:-1].split("\t")
        gene = d[0]
        nb_regions = int(d[1])
        isoform_ids = d[2].split(";")
        region_ids = d[3].split(";")
        regions = []
        for region in d[4:]:
            regions.append(ast.literal_eval(region))
        # Gene with all isoforms processed
        if nb_regions == len(region_ids):
            index_out.write(line)
            cat_regions = regions[:]
        # Gene with missing isoforms
        else:
            # Gene with non-matched isoforms
            if not scaffolds.has_key(gene):
                longest = ["", [], 0]
                # Get longest non-matched isoform
                if len(sequences[gene].keys()) != len(isoform_ids):
                    for iso in sorted(sequences[gene].keys()):
                        if not iso in isoform_ids:
                            tot_length = 0
                            ls_regions = []
                            for region in sorted(sequences[gene][iso]):
                                tot_length += len(sequences[gene][iso][region])
                                ls_regions.append(region)
                            if tot_length > longest[2]:
                                longest = [iso, ls_regions, tot_length]
                # Get longest non-matched region
                else:
                    for iso in sorted(sequences[gene].keys()):
                        for region in sorted(sequences[gene][iso]):
                            if not region in region_ids:
                                length = len(sequences[gene][iso][region])
                                if length > longest[2]:
                                    longest = [iso, [region], length]
                iso = longest[0]
                if not iso in isoform_ids:
                    isoform_ids.append(iso)
                for id in longest[1]:
                    region_ids.append(id)
                    regions.append([id, 1, len(sequences[gene][iso][id])])
                cat_regions = regions[:]
            # Gene with matched isoform
            else:
                cat_regions = []
                for region in regions:
                    test = True
                    s_isoform = region[0]
                    # Parsed region
                    if s_isoform in scaffolds[gene]:
                        for q_isoform in scaffolds[gene][s_isoform]:
                            isoform_id = q_isoform.split("|")[0]
                            # Modified region
                            if len(scaffolds[gene][s_isoform][q_isoform]) != 1:
                                test = False
                                cat_seq_id = s_isoform+"|"+q_isoform
                                s_seq = db[gene][s_isoform]
                                q_seq = sequences[gene][isoform_id][q_isoform]
                                cat_seq = ""
                                for s_block, q_block in scaffolds[gene][s_isoform][q_isoform]:
                                    # Block of the query sequence to include
                                    if s_block == []:
                                        q_start, q_stop = q_block[1:]
                                        cat_seq += q_seq[q_start-1:q_stop]
                                    # Block of the subject sequence to include
                                    else:
                                        s_start, s_stop = s_block[1:]
                                        cat_seq += s_seq[s_start-1:s_stop]
                            # Include isoform and region IDs
                            if not isoform_id in isoform_ids:
                                isoform_ids.append(isoform_id)
                            region_ids.append(q_isoform)
                    # Non-parsed or non-modified region
                    if test == True:
                        cat_regions.append(region)
                    # Update consensus region
                    else:
                        cat_regions.append([cat_seq_id, 1, len(cat_seq)])
                        db[gene][cat_seq_id] = cat_seq
            # Output updated index
            index_out.write("\t".join([gene, str(nb_regions), ";".join(sorted(isoform_ids)), ";".join(sorted(region_ids)), "\t".join(str(region) for region in cat_regions)])+"\n")
        # Output processed regions
        for region in cat_regions:
            region_id = region[0]
            isoform = region_id.split("|")[0]
            if db[gene].has_key(region_id):
                seq = db[gene][region_id]
            else:
                seq = sequences[gene][isoform][region_id]
            fasta_out.write(">"+gene+"|"+region_id+" "+"len="+str(len(seq))+"\n")
            for b in range((len(seq)/60)+1):
                text = seq[b*60:(b+1)*60]
                fasta_out.write(text)
                if text != "":
                    fasta_out.write("\n")
fasta_out.close()
index_out.close()
