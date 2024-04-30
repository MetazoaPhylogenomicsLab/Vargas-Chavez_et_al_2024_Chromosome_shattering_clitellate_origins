'''
3.3.recover_HOGs_from_lost.py
script to recoverer HOGs that are not lost in a minimum of species by node
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
    parser = argparse.ArgumentParser(description='Selecting duplicated and expanded hogs in an specific node')
    parser.add_argument('-q', type=str, help='path to directory with sequences lists')
    parser.add_argument('-d', type=str, help='list of LOST HOGs in the node (obtained from gene_rep_evol_at_hand.sh)')
    parser.add_argument('-l', type=str, help='list of species in the node')
    parser.add_argument('-p', type=float, help='cut-off for species number. float value') # use integer value if you want set the cutoff to a specific number of species
    parser.add_argument('-s', type=str, help='cut-off for sequences number. Integer value')
    parser.add_argument('-o', type=str, help='output file') #output file for recovered HOGs from LOST. THESE HOGS ARE NOT LOST!!!
    args = parser.parse_args()

    # Directory containing the sequences files
    directory = args.q

    #Set list of lost HOGs
    with open(args.d, 'r') as file:
    # Create an empty list to store the lines
        losses = []
        # Iterate over the lines of the file
        for line in file:
        # Remove the newline character at the end of the line
            lost = line.strip()
            # Append the line to the list
            losses.append(lost)

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
    spnumber = len(species)

    #Cut-off
    if args.p.is_integer():
        cutoff_sp = int(args.p)
    else:
        cutoff_sp = int(spnumber * float(args.p))
    print('this is the sp cutoff: ' + str(cutoff_sp))
    cutoff_seq = int(args.s)
    print('this is the seq cutoff: ' + str(cutoff_seq))
    
    for HOG in losses:
        #Obtain sequence file
        seqs_file = os.path.join(directory, HOG + '_seqs.txt')
        #Read the contents of the txt file
        with open(seqs_file, 'r') as file:
            elements = file.read().splitlines()
            elements = [item.split('.')[0] for item in elements]
            count_by_sp = Counter(elements)
            selected_elements = {element: count for element, count in count_by_sp.items() if count >= cutoff_seq and element in species}
            if len(selected_elements) >= cutoff_sp:
                with open(args.o, 'a') as ret:
                    ret.write(HOG + '\n')
                    ret.close()

if os.path.exists(args.o):
    with open(args.o, 'r') as file:
        HOGs_count = len(file.read().splitlines())
        print(HOGs_count, 'HOGs are not lost. These HOGs have at least', cutoff_seq, 'sequences in minimum', cutoff_sp, 'species of', spnumber, 'species included in the clade')
else:
    print('upsss!! there is no recovered lost genes in this node')
