#!/bin/bash

# 3.2.extract_filtered_GRE.sh
# script to obtain the filtered GRE from pyHam outputs
# wrote by Lisandra Benítez Álvarez on July 2023 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, Catalonia, Spain

RET=$1 #list of retained HOGs after filtering

for FOLDER in *_node_GRE
do
	cd $FOLDER
	for file in *_HOGs.txt
	do
		sed -i 's/HOG:/HOG/' $file
		grep -f $file ../$RET > ${file%%.txt}"_filtered.txt"
	done
	cd ..
done
