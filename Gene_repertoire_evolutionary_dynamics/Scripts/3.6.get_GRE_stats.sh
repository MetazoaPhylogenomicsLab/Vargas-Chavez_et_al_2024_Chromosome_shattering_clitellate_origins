#!/bin/bash

# 3.6.get_GRE_stats.sh
# script to obtain stats from GRE
# wrote by Lisandra Benítez Álvarez on February 2024 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, Catalonia, Spain

for FOLDER in *_node_GRE
do
    echo "************************" $FOLDER "*********************"
    cd $FOLDER
    wc -l GAINED_HOGs_filtered.txt | sed 's/^[ \t]*//' | awk '{print $2,$1}'
    wc -l REAL_LOST_filtered.txt | sed 's/^[ \t]*//' | awk '{print $2,$1}'
    wc -l REAL_RETAINED_HOGs_filtered.txt | sed 's/^[ \t]*//' | awk '{print $2,$1}'
    wc -l REAL_DUPLICATED_HOGs_filtered.txt | sed 's/^[ \t]*//' | awk '{print $2,$1}'
    wc -l REAL_EXPANDED_HOGs_filtered.txt | sed 's/^[ \t]*//' | awk '{print $2,$1}'
    echo "**********************************************************"
    cd ..
done
