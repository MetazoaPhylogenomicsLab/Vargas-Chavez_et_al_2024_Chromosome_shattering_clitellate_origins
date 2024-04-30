"""
reduce_HOG_annotation.py
script to reduce the annotation of HOGs
wrote by Lisandra Benítez Álvarez on August 2023 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, Catalonia, Spain
"""
import argparse
import time
import warnings
import pandas as pd
import numpy as np
import os
import shutil

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="select the GO terms by their frequency respect to the sequences number")
    parser.add_argument("-d", type=str, help="input folder with HOG annotation files by sequence", required=True)
    parser.add_argument("-p", type=float, help="Minimum percent of sequences where the GO has to be present", required=True)
    parser.add_argument("-o", type=str, help="path to out directory", required=True)
    parser.add_argument("-x", type=str, help="output file to save failed reduction", required=True)
    args = parser.parse_args()
    start_time = time.time()

    directory = args.d
    out_dir = args.o
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    file_names = [file for file in os.listdir(directory) if file.endswith('.txt')]
    
    for file_name in file_names:
        input_file = file_name
        output_file = input_file.replace("_anno_", "_HOG_anno_red_")
        
        HOG = input_file
        HOG = HOG.split('_anno_')[0]

        df = pd.read_table(f"{args.d}/{file_name}", header=None)
        df.set_index([0], inplace=True)

        sel = args.p
        #Here, we are setting the sel cut-off based on the value of args.p
        # args.p represent the cuttoff value to selec the GO by its frequency.
        #The GO frequency is calculated dividing the GO count by the number of sequences
        #If args.p = 0.5 we will select the GOs that are, at least in the 50% of the sequences
        #If args.p = 0.1 we will select the GOs that are, at least in the 10% of the sequences 

        df_countGO = df[1].value_counts()
        # print("df_countGO look like this")
        # print(df_countGO)
        countseqs = df.index.nunique()
        # print("Number of sequences")
        # print(countseqs)
        warnings.filterwarnings('ignore')
        df_sel = df_countGO/countseqs >= sel
        df_true = df_sel.loc[df_sel == True]
        # print("df_true look like this")
        # print(df_true)
        if df_true.size >= 1: 
            print(f"{HOG} has reduced ontology")
            df_true = pd.DataFrame(df_true)
            df_true = df_true.replace(True, HOG)
            df_true = df_true.reset_index()
            df_true.to_csv(f"{out_dir}/{output_file}.{args.p}", sep="\t", header = False, index = False)
        else:
            print(f"{input_file} does not have reduced ontology")
            with open(args.x, 'a+') as out_failed:
                out_failed.write(f"{input_file}.{args.p}\n")
    
    print(f"--- {time.time() - start_time} seconds ---")