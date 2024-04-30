#!/bin/sh

$IQTREE_PATH/iqtree -s Chapaev_frame_plus_RepeatPeps_both_virus_trimal_automated.aln -m MF -mset LG+F+G,WAG+F+G,JTT+F+G -madd EX2,EX3,EHO,LG4M,LG4X -redo

$IQTREE_PATH/iqtree -s Chapaev_frame_plus_RepeatPeps_both_virus_trimal_automated_shortseqs_removed.aln -m MF -mset LG+F+G,WAG+F+G,JTT+F+G -madd EX2,EX3,EHO,LG4M,LG4X -redo
