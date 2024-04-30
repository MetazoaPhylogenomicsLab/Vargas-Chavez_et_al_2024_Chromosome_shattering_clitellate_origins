#!/bin/sh

# 5.1.capture_genes_in_OMA_output.sh
# script to capture genes in OMA's output
# wrote by Lisandra Benítez Álvarez on February 2024 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, Catalonia, Spain

HOG=$1 #species composition table filtered HOGs
SP_HOGs=$2 #Position species in table HOGs

for file in *.txt
do

	SET=${file%%.txt}""
	echo "*************" analysing $SET "****************"

	countHOGs=$(grep -w -f $file $HOG | cut -f1,"$SP_HOGs" | wc -l)
	echo $countHOGs matching sequences agains HOGs
	grep -w -f $file $HOG | cut -f1,"$SP_HOGs" > $SET.match_HOGs.txt
	while read -r LINE;
	do
		matchHOG=$(grep -w "$LINE" $HOG | cut -f1)
		echo -e $LINE'\t'$matchHOG >> $SET.match_seqs_HOGs.txt
	done < $file

	echo "________________________________________________________________________________"

done

echo "*********************************************************************************"
echo "*******************************" SOME STATS "************************************"
echo "*********************************************************************************"

#Obtaining lists of identified HOGs
for file in *.match_HOGs.txt; do cut -f1 $file > ${file%%.txt}.list; done

#Obtaining lists of assigned and unassigned sequences
for file in *.match_seqs_HOGs.txt; do cat $file | awk 'BEGIN {FS="\t"} $2=="" {print}' > ${file%%.txt}.unassig.list; cat $file | awk 'BEGIN {FS="\t"} $2!="" {print}' > ${file%%.txt}.assig.list; done


#Obtaining some stats
echo Identified HOGs in UP
wc -l *UP*.match_HOGs.list | sed 's/^[ \t]*//' | grep ".list" | sed 's/.match_HOGs.list//' | awk '{print $2,$1}'
echo "*********************************************************************************"

echo Identified HOGs in Down
wc -l *Down*.match_HOGs.list | sed 's/^[ \t]*//' | grep ".list" | sed 's/.match_HOGs.list//' | awk '{print $2,$1}'
echo "*********************************************************************************"

echo Assigned sequences in UP
wc -l *UP*.match_seqs_HOGs.assig.list | sed 's/^[ \t]*//' | grep ".list" | sed 's/.match_seqs_HOGs.assig.list//' | awk '{print $2,$1}'
echo "*********************************************************************************"

echo Unassigned sequences in UP
wc -l *UP*.match_seqs_HOGs.unassig.list | sed 's/^[ \t]*//' | grep ".list" | sed 's/.match_seqs_HOGs.unassig.list//' | awk '{print $2,$1}'
echo "*********************************************************************************"

echo Assigned sequences in Down
wc -l *Down*.match_seqs_HOGs.assig.list | sed 's/^[ \t]*//' | grep ".list" | sed 's/.match_seqs_HOGs.assig.list//' | awk '{print $2,$1}'
echo "*********************************************************************************"

echo Unassigned sequences in Down
wc -l *Down*.match_seqs_HOGs.unassig.list | sed 's/^[ \t]*//' | grep ".list" | sed 's/.match_seqs_HOGs.unassig.list//' | awk '{print $2,$1}'
echo "*********************************************************************************"
