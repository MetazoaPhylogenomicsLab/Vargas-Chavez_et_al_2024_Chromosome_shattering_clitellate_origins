'''
2.3.filtering_HOGs_by_CompletenessScore.py
script to filter HOGs from OMA's output
wrote by Lisandra Benítez Álvarez on July 2023 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, Catalonia, Spain
'''
import os
import pandas as pd
from collections import Counter
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filtering HOGs by Completeness Score')
    parser.add_argument('-d', type=str, help='path to directory with sequences lists')
    parser.add_argument('-l', type=str, help='list of species by node')
    parser.add_argument('-c', type=str, help='cut-off for filtering as float Ex. 0.2')
    args = parser.parse_args()

    # Directory containing the sequences files
    directory = args.d

    #List of species by node
    list = args.l
    #Import list as dataframe
    df_list = pd.read_csv(list, delimiter= '\t', index_col = "node")

    #Cut-off
    cutoff = args.c
    cutoff = float(cutoff)

    # Get a list of txt files in the directory
    seqs_files = [filename for filename in os.listdir(directory) if filename.endswith('.txt')]

    # Loop through each txt file
    for filename in seqs_files:
        file_path = os.path.join(directory, filename)
        HOG = filename.replace("_seqs.txt", "")
        node = filename.split('_')[1]
        node = int(node)
        sp_node = df_list.loc[node]['species']
        sp_node = int(sp_node)
        score = sp_node * cutoff
        score = int(score)

        # Read the contents of the txt file
        with open(file_path, 'r') as file:
            elements = file.read().splitlines()
            elements = [item.split('_')[0] for item in elements]
            count_by_sp = Counter(elements)
            total_count = len(count_by_sp)
            if total_count > score:
                print(HOG, 'has passed the filtering.', HOG, 'has been saved in retained_HOGs.txt')
                with open('retained_HOGs.txt', 'a') as ret:
                    ret.write(HOG + '\n')
                    ret.close()
            else:
                print(HOG, 'has not passed the filtering.', HOG, 'has been saved in discarted_HOGs.txt')
                with open('discarted_HOGs.txt', 'a') as dis:
                    dis.write(HOG + '\n')
                    dis.close()

with open('retained_HOGs.txt', 'r') as file:
    HOGs_count = len(file.read().splitlines())
    print(HOGs_count, 'HOGs have been retained after filtering')

with open('discarted_HOGs.txt', 'r') as file:
    HOGs_count = len(file.read().splitlines())
    print(HOGs_count, 'HOGs have been discarted after filtering')


