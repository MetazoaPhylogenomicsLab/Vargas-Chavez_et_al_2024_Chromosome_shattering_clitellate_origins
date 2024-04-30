"""
4.1.annotate_OMA_HOGs_with_gopredsim.py
script to transfer the annotation of sequences done with Fantasia to HOGs inferred with OMA 
wrote by Lisandra Benítez Álvarez on August 2023 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, Catalonia, Spain
"""
import argparse
import time
import textwrap
import os
import pandas as pd

ANNOTATIONS_TYPES = ["bpo", "cco", "mfo"]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate HOGs using the sequence's annotation done with gopredsim")
    parser.add_argument("-s", type=str, help="input folder with HOGs sequence lists", required=True)
    parser.add_argument("-a", type=str, help="input folder with species annotations", required=True)
    parser.add_argument("-l", type=str, help="species list", required=True)
    parser.add_argument("-o", type=str, help="path to out directory", required=True)
    args = parser.parse_args()
    start_time = time.time()
    
    dir_seq = args.s
    dir_anno = args.a
    out_dir = args.o
    os.makedirs(out_dir)

    with open(args.l, "r") as file:
        sp_list = [line.strip() for line in file]

    annotations_df = {ANNOTATIONS_TYPES[0]: None, ANNOTATIONS_TYPES[1]: None, ANNOTATIONS_TYPES[2]: None}
    for species in sp_list:
        for anno_type in ANNOTATIONS_TYPES:
            anno_df = pd.read_table(
                f"{args.a}/gopredsim_{species}_prott5_1_{anno_type}.txt", #f"{args.a}/{species}_anno_{anno_type}.txt",
                header = None,
                index_col=0,
            )
            annotations_df[anno_type] = pd.concat([annotations_df[anno_type], anno_df])
   
    file_names = [file for file in os.listdir(dir_seq) if file.endswith('_seqs.txt')]
    
    for file_name in file_names:
        seq_file = file_name
        HOG = seq_file.replace("_seqs.txt", "")
        with open(f"{dir_seq}/{seq_file}", "r") as f:
            HOG_seqs = [line.strip() for line in f]
            #print(HOG_seqs)
        df_seq_anno_bpo = []
        for seq in HOG_seqs:
            try:
                seq_anno_bpo = annotations_df[ANNOTATIONS_TYPES[0]].loc[[seq]]
                #print(seq_anno_bpo)
            except KeyError:
                #print(f"{seq} without {ANNOTATIONS_TYPES[0]} annotation, the annotation of an identical sequence will be propagated")
                continue
            df_seq_anno_bpo.append(seq_anno_bpo)
        if len(df_seq_anno_bpo) > 0:
            df_seq_anno_bpo = pd.concat(df_seq_anno_bpo)
            #print(df_seq_anno_bpo)
            df_seq_anno_bpo.to_csv(
            f"{args.o}/{HOG}_anno_{ANNOTATIONS_TYPES[0]}.txt",
            sep="\t",
            index_label=False,
            header=False,
            )
        df_seq_anno_cco = []
        for seq in HOG_seqs:
            try:
                seq_anno_cco = annotations_df[ANNOTATIONS_TYPES[1]].loc[[seq]]
            except KeyError:
                #print(f"{seq} without {ANNOTATIONS_TYPES[1]} annotation, the annotation of an identical sequence will be propagated")
                continue
            df_seq_anno_cco.append(seq_anno_cco)
        if len(df_seq_anno_cco) > 0:
            df_seq_anno_cco = pd.concat(df_seq_anno_cco)
            df_seq_anno_cco.to_csv(
            f"{args.o}/{HOG}_anno_{ANNOTATIONS_TYPES[1]}.txt",
            sep="\t",
            index_label=False,
            header=False,
            )
        df_seq_anno_mfo = []
        for seq in HOG_seqs:
            try:
                seq_anno_mfo = annotations_df[ANNOTATIONS_TYPES[2]].loc[[seq]]
            except KeyError:
                #print(f"{seq} without {ANNOTATIONS_TYPES[2]} annotation, the annotation of an identical sequence will be propagated")
                continue
            df_seq_anno_mfo.append(seq_anno_mfo)
        if len(df_seq_anno_mfo) > 0:
            df_seq_anno_mfo = pd.concat(df_seq_anno_mfo)
            df_seq_anno_mfo.to_csv(
            f"{args.o}/{HOG}_anno_{ANNOTATIONS_TYPES[2]}.txt",
            sep="\t",
            index_label=False,
            header=False,
            )

print(f"Process time elapsed: {time.time() - start_time} seconds.")