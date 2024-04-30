#!/bin/bash

# 3.3.launch_recover_HOGs_from_lost.sh
# script to launch 3.3.recover_HOGs_from_lost.py script in several folders
# wrote by Lisandra Benítez Álvarez on February 2024 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, Catalonia, Spain

SCRIPT=$1 #path to script 3.3.recover_HOGs_from_lost.py
HOGS_SEQs=$2 #path to folder with sequences by HOG

for FOLDER in *_node_GRE
do
    echo "************************" $FOLDER "*********************"
    cd $FOLDER

    python3 ../$SCRIPT -q ../$HOGS_SEQs -d=LOST_par_to_anc_HOGs_filtered.txt -l=species_by_node.txt -p 0.2 -s 1 -o=recovered_from_LOST_par_to_anc_HOGs_filtered.txt
    grep -v -f recovered_from_LOST_par_to_anc_HOGs_filtered.txt LOST_par_to_anc_HOGs_filtered.txt > REAL_LOST_filtered.txt
    cat RETAINED_HOGs_filtered.txt recovered_from_LOST_par_to_anc_HOGs_filtered.txt >> REAL_RETAINED_HOGs_filtered.txt

    echo HOGs recovered from previous LOST
    wc -l recovered_from_LOST_par_to_anc_HOGs_filtered.txt

    echo new LOST after recover
    wc -l REAL_LOST_filtered.txt

    echo new RETAINED after recover
    wc -l REAL_RETAINED_HOGs_filtered.txt

    cd ..
done
