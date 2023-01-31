#!/usr/bin/env python

# Bi623 - Human and Zebrafish Reciprocal Best Hits
# Using longest protein data from Bi621 Problem Set 7
# VARIABLE NAMES: h/human or z/zfish indicate the QUERY, database is the other 

def mart_dict(fname: str):
    '''
    Takes a path to a tab-separated Biomart export file containing columns:
        Gene stable ID	Protein stable ID	Gene name
    Returns a dictionary with:
        keys = protein stable ID, values = (gene stable ID, gene name)
    Excludes records with no value for Protein Stable ID.
    '''
    lookup: dict = {}
    names: list = []

    with open(fname, 'r') as file:
        for line in file:
            names = line.strip('\n').split('\t')
            lookup[names[1]] = (names[0], names[2])
    # Remove header line and blank protein entry before returning.
    lookup.pop("Protein stable ID")
    lookup.pop("")

    return lookup

def best_hit_filter(sorted_path: str) -> dict:
    '''
    Takes a path to a BLASTp .tsv output file.
    File must be sorted by qseqid and then by evalue!
    Returns a dictionary of best hits for each protein query:
        keys = prot query ID, values = [match prot ID, evalue, discard flag]
    Discard flag is TRUE if the protein query had multiple best hits.
    '''
    best_hits: dict = {}
    with open(sorted_path, 'r') as fname:
        for line in fname:
            columns: list = line.strip().split()
            protid: str = columns[0]
            matchid: str = columns[1]
            evalue: float = float(columns[10])

            if protid in best_hits:
                # If new match has lower evalue, save the new match, reset discard flag.
                if evalue < best_hits[protid][1]:
                    best_hits[protid] = [matchid, evalue, False]

                # If new match has the same evalue as saved evalue, mark entry for discard.
                elif evalue == best_hits[protid][1]:
                    best_hits[protid][2] = True

            # If hit is not already saved, save it.
            else:
                best_hits[protid] = [matchid, evalue, False]

    # Remove all entries marked for discard.
    for key in list(best_hits.keys()):
        if best_hits[key][2] == True:
            best_hits.pop(key)

    return best_hits


# Create lookup tables from biomart export files.
h_mart: str = "/projects/bgmp/bmeluch/bioinfo/Bi621/PS/ps7-bmeluch/biomart_export/human_mart_expt.txt"
z_mart: str = "/projects/bgmp/bmeluch/bioinfo/Bi621/PS/ps7-bmeluch/biomart_export/zebrafish_mart_expt.txt"

zfish_lookup = mart_dict(z_mart)
human_lookup = mart_dict(h_mart)

# Create dictionaries of best hits from BLASTp output files.
h_sort_path = "/projects/bgmp/bmeluch/bioinfo/Bi623/PS/PS1/blastp_output_sorted/h_to_zdb_sort.tsv"
z_sort_path = "/projects/bgmp/bmeluch/bioinfo/Bi623/PS/PS1/blastp_output_sorted/z_to_hdb_sort.tsv"

h_best_hits: dict = best_hit_filter(h_sort_path)
z_best_hits: dict = best_hit_filter(z_sort_path)


# Create empty dictionary for reciprocal best hits.
# Key = human protID, value = zfish protID
RBH: dict = {}

# Iterate through zfish best hits dictionary.
for zprot in z_best_hits:
    # Check if the human match for zprot query is a key in the human dict.
    if z_best_hits[zprot][0] in h_best_hits:
        # Check reciprocity, save match if the proteins are the same.
        if h_best_hits[z_best_hits[zprot][0]][0] == zprot:
            RBH[z_best_hits[zprot][0]] = zprot
print("RBH: ", len(RBH))

with open("/projects/bgmp/bmeluch/bioinfo/Bi623/PS/PS1/Human_Zebrafish_RBH.tsv", 'w') as tsv:
    tsv.write("Human Gene ID\tHuman Protein ID\tHuman Gene Name\t"+
    "Zebrafish Gene ID\tZebrafish Protein ID\tZebrafish Gene Name\n")
    # Use protein codes to look up gene information from biomart export.
    for k in RBH:
        tsv.write(human_lookup[k][0]+'\t'+
        k+'\t'+
        human_lookup[k][1]+'\t'+
        zfish_lookup[RBH[k]][0]+'\t'+
        RBH[k]+'\t'+
        zfish_lookup[RBH[k]][1]+'\n')
