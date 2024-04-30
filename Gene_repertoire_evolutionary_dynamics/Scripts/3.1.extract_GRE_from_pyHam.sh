#!/bin/sh

# 3.1.extract_GRE_from_pyHam.sh
# script to obtain the first aproximation of the GRE from pyHam outputs
# wrote by Lisandra Benítez Álvarez on July 2023 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, Catalonia, Spain

NODE_LIST=$1 # File to run GRE automatically. See file in EXAMPLE folder
TOP_HOG=$2 #Top hog list top_hog_pyham.txt
ANC_DIR=$3 #Directory containing ancestral genomes from pyHam

while read -r LINE
do
    NODE=$(echo $LINE | cut -d " " -f1) #Node to analyze
    MY_GENOME=$(echo $LINE | cut -d " " -f2)        #Ancestral genome of interest obtained from pyHam 
    PARENTAL_GENOME=$(echo $LINE | cut -d " " -f3)  #Parental genome obtained from pyHam
    TAXON=$(echo $LINE | cut -d " " -f4) #Taxon id of interest node         
    TAXONP=$(echo $LINE | cut -d " " -f5) #Taxon id parent node
    OUT_DIR=${NODE}_node_GRE

    mkdir $OUT_DIR

    echo "*************************************************"  > $OUT_DIR/report.txt
    echo "********** GENE REPERTORIE EVOLUTION ************" >> $OUT_DIR/report.txt
    echo "*************************************************" >> $OUT_DIR/report.txt 

    echo "Total length of the ancestral genome"              >> $OUT_DIR/report.txt
    T_len=$(cat $ANC_DIR/$MY_GENOME | tr ',' '\n' | wc -l)
    echo $T_len                                              >> $OUT_DIR/report.txt
    echo "*************************************************" >> $OUT_DIR/report.txt

    echo "Identified HOGs in the ancestral genome"           >> $OUT_DIR/report.txt
    LEN=$(cat $ANC_DIR/$MY_GENOME | tr ',' '\n' | grep -v "<HOG()>" | wc -l) 
    echo $LEN                                                >> $OUT_DIR/report.txt
    echo "*************************************************" >> $OUT_DIR/report.txt

    echo "Gained HOGs in the ancestral genome"               >> $OUT_DIR/report.txt
    GAIN=$(cat $TOP_HOG | tr ',' '\n' | grep "_$TAXON" | cut -d " " -f3 | sed 's/<HOG(//;s/)>//' | wc -l)
    echo $GAIN                                               >> $OUT_DIR/report.txt 
    cat $TOP_HOG | tr ',' '\n' | grep "_$TAXON" | cut -d " " -f3 | sed "s/<HOG(//;s/_$TAXON)>//" > $OUT_DIR/GAINED_HOGs.txt
    echo "*************************************************" >> $OUT_DIR/report.txt

    echo "Total duplications events in the ancestral genome" >> $OUT_DIR/report.txt
    DUP_ev=$(cat $ANC_DIR/$MY_GENOME | grep -oE "<HOG\(HOG:[0-9]+(\.[0-9a-z]+)+_$TAXON)>" | sed 's/<HOG(//;s/)>//' | wc -l)
    echo $DUP_ev                                             >> $OUT_DIR/report.txt
    echo "Total duplicated HOGs in the ancestral genome"     >> $OUT_DIR/report.txt
    DUP=$(cat $ANC_DIR/$MY_GENOME | grep -oE "<HOG\(HOG:[0-9]+(\.[0-9a-z]+)+_$TAXON)>" | sed 's/<HOG(//;s/)>//' | cut -d "." -f1 | sort | uniq | wc -l)
    echo $DUP                                                >> $OUT_DIR/report.txt
    cat $ANC_DIR/$MY_GENOME | grep -oE "<HOG\(HOG:[0-9]+(\.[0-9a-z]+)+_$TAXON)>" | sed 's/<HOG(//;s/)>//' | cut -d "." -f1 | sort | uniq > $OUT_DIR/DUPLICATED_HOGs.txt
    echo "*************************************************" >> $OUT_DIR/report.txt

    echo "Retained HOGs in the ancestral genome"             >> $OUT_DIR/report.txt
    RET=$(cat $ANC_DIR/$MY_GENOME | tr ',' '\n' | grep -v "<HOG()>" | grep -v -f $OUT_DIR/GAINED_HOGs.txt | sed 's/<HOG(//;s/)>//' | sed "s/_$TAXON//" | cut -d "." -f1 | sort | uniq | wc -l)
    echo $RET                                                >> $OUT_DIR/report.txt
    cat $ANC_DIR/$MY_GENOME | tr ',' '\n' | grep -v "<HOG()>" | grep -v -f $OUT_DIR/GAINED_HOGs.txt | sed 's/<HOG(//;s/)>//' | sed "s/_$TAXON//" | cut -d "." -f1 | sort | uniq | sed 's/^ *//' > $OUT_DIR/RETAINED_HOGs.txt
    echo "*************************************************" >> $OUT_DIR/report.txt

    echo "Calculating gene lost from parent node"            >> $OUT_DIR/report.txt
    echo "Total length of the parent genome"                 >> $OUT_DIR/report.txt
    T_lenP=$(cat $ANC_DIR/$PARENTAL_GENOME | tr ',' '\n' | wc -l)
    echo $T_lenP                                             >> $OUT_DIR/report.txt
    echo "Identified HOGs in the parent genome"              >> $OUT_DIR/report.txt
    LENP=$(cat $ANC_DIR/$PARENTAL_GENOME | tr ',' '\n' | grep -v "<HOG()>" | wc -l)
    echo $LENP                                               >> $OUT_DIR/report.txt
    echo "*************************************************" >> $OUT_DIR/report.txt

    echo "Lost genes in the ancestral genome respect to the parent node" >> $OUT_DIR/report.txt
    echo "Lost HOGs"                                                     >> $OUT_DIR/report.txt
    LOST=$(cat $ANC_DIR/$PARENTAL_GENOME | tr ',' '\n' | grep -v "<HOG()>" | grep -v -f $OUT_DIR/RETAINED_HOGs.txt | sed "s/_$TAXONP)>//" | cut -d "." -f 1 | sort | uniq | cut -d "(" -f2 | wc -l)
    echo $LOST                                                           >> $OUT_DIR/report.txt
    cat $ANC_DIR/$PARENTAL_GENOME | tr ',' '\n' | grep -v "<HOG()>" | grep -v -f $OUT_DIR/RETAINED_HOGs.txt | sed "s/_$TAXONP)>//" | cut -d "." -f 1 | sort | uniq | cut -d "(" -f2 > $OUT_DIR/LOST_par_to_anc_HOGs.txt

done < $NODE_LIST
