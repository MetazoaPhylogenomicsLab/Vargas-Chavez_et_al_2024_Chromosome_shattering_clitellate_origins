#!/bin/bash

# 3.8.extract_seqs_from_GRE.sh
# script to extract the sequences included in the HOGs of every set of the GRE
# wrote by Lisandra Benítez Álvarez on February 2024 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, Catalonia, Spain

SEQS=$1 # folder containing the seqs by HOG files

for FOLDER in node*
do
    cd $FOLDER
    
    for file in *_HOG.txt
    do
        while read -r LINE
        do
            HOG=$LINE
            cat $SEQS/${HOG}_seqs.txt >> ${file%%_HOG.txt}"_seqs.txt"
        done < $file 
    done

    cd ..
done
