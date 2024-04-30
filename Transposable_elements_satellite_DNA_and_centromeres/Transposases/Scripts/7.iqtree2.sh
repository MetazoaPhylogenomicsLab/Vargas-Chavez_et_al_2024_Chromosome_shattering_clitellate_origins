#!/bin/sh

MODEL="WAG+F+R5"
$IQTREE_PATH/iqtree2 -s Chapaev_frame_plus_RepeatPeps_both_virus_trimal_automated.aln -m $MODEL -B 1000 -T 1 -alrt 1000 -redo

MODEL="WAG+F+R5"
$IQTREE_PATH/iqtree2 -s Chapaev_frame_plus_RepeatPeps_both_virus_trimal_automated_shortseqs_removed.aln -m $MODEL -B 1000 -T 1 -alrt 1000 -redo
