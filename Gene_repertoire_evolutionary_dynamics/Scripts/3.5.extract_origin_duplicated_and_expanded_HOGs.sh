#!/bin/bash

# 3.5.extract_origin_duplicated_and_expanded_HOGs.sh
# script to extract the origin of duplicated and expanded HOGs
# wrote by Lisandra Benítez Álvarez on February 2024 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, Catalonia, Spain

for FOLDER in *_node_GRE
do
	echo "************************" obtaining origin for duplicated and expanded HOGs in node ${FOLDER%%_node_GRE}"" "**********************"
    cd $FOLDER
        
		echo origin for duplicated HOGs
		for file in REAL_DUPLICATED_HOGs_filtered.txt*
		do
			cat $file | cut -d "_" -f2 | sort | uniq -c | sed 's/^[ \t]*//' | awk '{print $2,$1}' > ${file%%.txt}"_origin.txt"
			cat ${file%%.txt}"_origin.txt"
		done
        
		echo origin for expanded HOGs
		for file in REAL_EXPANDED_HOGs_filtered.txt*
		do
			cat $file | cut -d "_" -f2 | sort | uniq -c | sed 's/^[ \t]*//' | awk '{print $2,$1}' > ${file%%.txt}"_origin.txt"
			cat ${file%%.txt}"_origin.txt"
		done
	
	cd ..
done
