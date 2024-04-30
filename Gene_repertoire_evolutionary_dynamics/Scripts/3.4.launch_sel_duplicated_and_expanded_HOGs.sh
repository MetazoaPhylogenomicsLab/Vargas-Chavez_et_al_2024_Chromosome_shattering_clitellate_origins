#!/bin/bash

# 3.4.launch_sel_duplicated_and_expanded_HOGs.sh
# script to launch 3.4.sel_duplicated_and_expanded_HOGs.py script
# wrote by Lisandra Benítez Álvarez on February 2024 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, Catalonia, Spain

SCRIPT=$1 #path to script 3.4.sel_duplicated_and_expanded_HOGs.py
HOGS_SEQs=$2 #path to folder with sequences by HOG
PDUP=$3 #recommendable 0.2
SDUP=$4 #recommendable 2
PEXP=$5 #recommendable 0.5
SEXP=$6 #recommendable 3

for FOLDER in *_node_GRE
do
    echo "************************" $FOLDER "*********************"
    cd $FOLDER
    
    python3 $SCRIPT -q ../$HOGS_SEQs -d=DUPLICATED_HOGs_filtered.txt -l=species_by_node.txt -p $PDUP -s $SDUP -o=REAL_DUPLICATED_HOGs_filtered.txt
    echo HOGs REAL DUPLICATED
    wc -l REAL_DUPLICATED_HOGs_filtered.txt

    python3 $SCRIPT -q ../$HOGS_SEQs -d=DUPLICATED_HOGs_filtered.txt -l=species_by_node.txt -p $PEXP -s $SEXP -o=REAL_EXPANDED_HOGs_filtered.txt
    echo HOGs REAL EXPANDED
    wc -l REAL_EXPANDED_HOGs_filtered.txt

    cd ..
done
