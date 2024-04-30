#!/bin/bash

# 3.7.copy_GRE2folder.sh
# script to copy final GRE (filtered and curated) to independent folders by node 
# wrote by Lisandra Benítez Álvarez on February 2024 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, Catalonia, Spain

OUT=$1 #folder to copy the files

for FOLDER in *_node_GRE
do
    NODE=${FOLDER%%_node_GRE}""
    mkdir $OUT/node${NODE}
    
    cd $FOLDER

    cp GAINED_HOGs_filtered.txt $OUT/node${NODE}/GAINED_node${NODE}_HOG.txt
    cp REAL_DUPLICATED_HOGs_filtered.txt $OUT/node${NODE}/DUPLICATED_node${NODE}_HOG.txt
    cp REAL_EXPANDED_HOGs_filtered.txt $OUT/node${NODE}/EXPANDED_node${NODE}_HOG.txt
    cp REAL_LOST_filtered.txt $OUT/node${NODE}/LOST_node${NODE}_HOG.txt
    cp REAL_RETAINED_HOGs_filtered.txt $OUT/node${NODE}/RETAINED_node${NODE}_HOG.txt
    
    cd ..
done        

