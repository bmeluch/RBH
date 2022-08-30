#!/usr/bin/env python

# VARIABLE NAMES: h/human or z/zfish indicate the QUERY not the database. Database is the other organism. 

#################### FUNCTIONS ####################

# save mart files as dictionaries
# key: protein stable id; value: (gene stable id, gene name)
def mart_dict(fname: str):
    '''Takes a tab-separated biomart export file containing columns 
    Gene stable ID	Protein stable ID	Gene name
    Returns a dictionary with 
    keys = protein stable ID, values = (gene stable ID, gene name)
    Excludes records with no value for Protein Stable ID.
    '''
    lookup: dict = {}
    names: list = []

    with open(fname, 'r') as file:
        for line in file:
            names = line.strip('\n').split('\t')
            lookup[names[1]] = (names[0], names[2])
    # cleanup: delete header line, blank protein entry
    lookup.pop("Protein stable ID")
    lookup.pop("")

    return lookup


# read in BLASTp .tsv output and only retain the best hits for each protein query based on e-value
# there may be multiple best hits with the same e-value. in that case, do not retain record for that protein query
def best_hit_filter(sorted_path: str) -> dict:
    '''Takes a path to a BLASTp .tsv output file containing standard columns 
    'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    File must be sorted by qseqid and then by evalue!
    Returns a dictionary with 
    keys = protein query ID, values = [match protein ID, evalue, discard flag]
    Discard flag is TRUE if the protein query had multiple best hits and should be disregarded.
    '''
    best_hits: dict = {}
    with open(sorted_path, 'r') as fname:
        for line in fname:
            columns: list = line.strip().split()
            protid: str = columns[0]
            matchid: str = columns[1]
            evalue: float = float(columns[10])

            # check if hit is already in dictionary
            if protid in best_hits:

                # if new match is better (lower evalue), record the new match
                # change discard flag back to false if a better match is found
                if evalue < best_hits[protid][1]:
                    best_hits[protid] = [matchid, evalue, False]

                # if new match has the same evalue as saved evalue, mark this entry to be discarded
                elif evalue == best_hits[protid][1]:
                    best_hits[protid][2] = True

            # if hit is not already saved, save it
            else:
                best_hits[protid] = [matchid, evalue, False]

    # remove all discard true entries
    for key in list(best_hits.keys()):
        if best_hits[key][2] == True:
            best_hits.pop(key)

    return best_hits


#################### READ IN FILES ####################

# import biomart output for lookup
h_mart: str = "/projects/bgmp/bmeluch/bioinfo/Bi621/PS/ps7-bmeluch/biomart_export/human_mart_expt.txt"
z_mart: str = "/projects/bgmp/bmeluch/bioinfo/Bi621/PS/ps7-bmeluch/biomart_export/zebrafish_mart_expt.txt"

zfish_lookup = mart_dict(z_mart)
human_lookup = mart_dict(h_mart)

# Sorted BLASTp output file paths
#'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
h_sort_path = "/projects/bgmp/bmeluch/bioinfo/Bi623/PS/PS1/blastp_output_sorted/h_to_zdb_sort.tsv"
z_sort_path = "/projects/bgmp/bmeluch/bioinfo/Bi623/PS/PS1/blastp_output_sorted/z_to_hdb_sort.tsv"

# unsorted BLASTp output file paths
# got it to work on unsorted files, take that, computer
# h_sort_path = "/projects/bgmp/bmeluch/bioinfo/Bi621/PS/ps7-bmeluch/results/humanp_to_zfishdb.tsv"
# z_sort_path = "/projects/bgmp/bmeluch/bioinfo/Bi621/PS/ps7-bmeluch/results/zfishp_to_humandb.tsv"

# test file paths
# h_sort_path = "/projects/bgmp/bmeluch/bioinfo/Bi623/PS/PS1/test_files/h_to_zdb_sort_test.tsv"
# z_sort_path = "/projects/bgmp/bmeluch/bioinfo/Bi623/PS/PS1/test_files/z_to_hdb_sort_test.tsv"

# dicts: keys = protein ID, values = [match protein ID, evalue, discard flag]
h_best_hits: dict = best_hit_filter(h_sort_path)
z_best_hits: dict = best_hit_filter(z_sort_path)

# check
print("Human to zfish database best hits: ", len(h_best_hits))
print("Zfish to human database best hits: ", len(z_best_hits))

#################### RECIPROCAL BEST HITS ####################

# key = human code, values = zebrafish code
RBH: dict = {}

# iterate through zfish best hits dictionary
for zprot in z_best_hits:
    # check if the human match for zprot query is in the human dictionary as a key(query)
    if z_best_hits[zprot][0] in h_best_hits:
        # if so, check if the zebrafish match for the human protein is the same as zprot
        if h_best_hits[z_best_hits[zprot][0]][0] == zprot:
            # if so, save the match to RBH
            RBH[z_best_hits[zprot][0]] = zprot
print("RBH: ", len(RBH))

# write out best hits
with open("/projects/bgmp/bmeluch/bioinfo/Bi623/PS/PS1/Human_Zebrafish_RBH.tsv", 'w') as tsv:
    tsv.write("Human Gene ID"+'\t'+"Human Protein ID"+'\t'+"Human Gene Name"+'\t'+
    "Zebrafish Gene ID"+'\t'+"Zebrafish Protein ID"+'\t'+"Zebrafish Gene Name"+'\n')
    # use rbh dictionary entries to look up gene information from biomart export
    for k in RBH:
        tsv.write(human_lookup[k][0]+'\t'+
        k+'\t'+
        human_lookup[k][1]+'\t'+
        zfish_lookup[RBH[k]][0]+'\t'+
        RBH[k]+'\t'+
        zfish_lookup[RBH[k]][1]+'\n')