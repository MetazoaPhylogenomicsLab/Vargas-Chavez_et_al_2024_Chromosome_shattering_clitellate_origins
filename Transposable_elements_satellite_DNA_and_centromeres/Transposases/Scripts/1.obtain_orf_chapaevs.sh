#!/bin/bash

#Input files and paths
WD=$STORE2/Annelida_project/transposases
CHAPAEV_FILE=$WD/CMC-Chapaev-3-families.fa
CHAPAEV_DB=$WD/Chapaev_RepeatPeps.fasta

#1) Build Chapaev prots database
makeblastdb -in $CHAPAEV_DB -dbtype prot

#2) Execute blastx using our transposases DNA fasta as query and the database as target
blastx -query $CHAPAEV_FILE -db $CHAPAEV_DB -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe" -out $WD/CMC-Chapaev-3-families_frame.blastout
