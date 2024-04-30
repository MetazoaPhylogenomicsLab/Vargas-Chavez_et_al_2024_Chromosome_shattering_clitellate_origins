'''
3.4.sel_duplicated_and_expanded_HOGs.py
script to obtain duplicated and expanded HOGs
wrote by Lisandra Benítez Álvarez on February 2024 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, Catalonia, Spain
'''
from itertools import count
import os
import pandas as pd
from collections import Counter
import argparse
import re
import glob

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Selecting duplicated and expanded HOGs in an specific node')
    parser.add_argument('-q', type=str, help='path to directory with sequences lists')
    parser.add_argument('-d', type=str, help='list of all DUPLICATED HOGs in the node (inferred from pyHam in the first approach)')
    parser.add_argument('-l', type=str, help='list of species in the node')
    parser.add_argument('-p', type=float, help='cut-off for species number. float value')
    parser.add_argument('-s', type=str, help='cut-off for sequences number. Integer value')
    parser.add_argument('-o', type=str, help='output file')
    args = parser.parse_args()

    # Directory containing the sequences files
    directory = args.q

    #Set list of duplicated HOGs
    with open(args.d, 'r') as file:
    # Create an empty list to store the lines
        duplications = []
        # Iterate over the lines of the file
        for line in file:
        # Remove the newline character at the end of the line
            duplicated = line.strip()
            # Append the line to the list
            duplications.append(duplicated)
    #print(duplications)

    #Set list of species
    with open(args.l, 'r') as file:
    # Create an empty list to store the lines
        species = []
        # Iterate over the lines of the file
        for line in file:
        # Remove the newline character at the end of the line
            specie = line.strip()
            # Append the line to the list
            species.append(specie)
    #print(species)
    spnumber = len(species)
    print('species number in node: ' + str(spnumber))

    #Cut-off
    if args.p.is_integer():
        cutoff_sp = int(args.p)
    else:
        cutoff_sp = int(spnumber * float(args.p))
    print('this is the sp cutoff: ' + str(cutoff_sp))
    cutoff_seq = int(args.s)

    for HOG in duplications:
        seqs_file = os.path.join(directory, HOG + '_seqs.txt')
        #Read the contents of the txt file
        with open(seqs_file, 'r') as file:
            elements = file.read().splitlines()
            elements = [item.split('.')[0] for item in elements] # use '_' for Platy and '.' for Annelida
            count_by_sp = Counter(elements)
            selected_elements = {element: count for element, count in count_by_sp.items() if count >= cutoff_seq and element in species}
            if len(selected_elements) >= cutoff_sp:
                with open(args.o, 'a') as ret:
                    ret.write(HOG + '\n')
                    ret.close()

if os.path.exists(args.o):
    with open(args.o, 'r') as file:
        HOGs_count = len(file.read().splitlines())
        print(HOGs_count, 'HOGs have been duplicated. These HOGs have at least', cutoff_seq, 'sequences in minimum', cutoff_sp, 'species of', spnumber, 'species included in the clade')
else:
    print('upsss!! there is no "real" duplicated genes in this node')
