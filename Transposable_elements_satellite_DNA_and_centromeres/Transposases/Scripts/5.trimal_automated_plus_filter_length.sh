#!/bin/sh
#Trimm alignment
$TRIMAL_PATH/trimal -keepheader -in Chapaev_frame_plus_RepeatPeps_both_virus.aln -automated1 > Chapaev_frame_plus_RepeatPeps_both_virus_trimal_automated.aln

cat Chapaev_frame_plus_RepeatPeps_both_virus_trimal_automated.aln | seqkit seq --remove-gaps -m 100 | grep ">" | tr -d ">" > sequences_larger100_after_trimal.txt

seqtk subseq Chapaev_frame_plus_RepeatPeps_both_virus_trimal_automated.aln sequences_larger100_after_trimal.txt > Chapaev_frame_plus_RepeatPeps_both_virus_trimal_automated_shortseqs_removed.aln
