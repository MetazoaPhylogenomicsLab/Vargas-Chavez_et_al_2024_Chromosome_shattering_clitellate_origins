# Transcriptomes

## 1. Description 
The transcriptome for *Hirudo medicinalis* after being subjected to several abiotic stress experiments: visible and UV-B light, hyperoxia, hypoxia, and osmoregulation, plus control specimens. All samples were dissected into head, mid body and tail, flash frozen in liquid nitrogen and kept <-70ºC until further processing. RNA extractions were performed using the TRIzol® reagent (Invitrogen, USA) method following the manufacturer’s instructions. Concentration of all samples was assessed by Qubit RNA BR Assay kit (Thermo Fisher Scientific). Samples were subjected to either Illumina’s TruSeq Stranded mRNA library preparation kit or Truseq stranded Total RNA Library with Ribo-Zero depending on the tissue and species, and sequenced on a NovaSeq 6000 (Illumina, 2 × 150 bp) for a 6Gb coverage. Quality assessment of the raw Illumina RNA-seq reads was performed using FastQC v0.11.9 (https://www.bioinformatics.babraham.ac.uk/projects/fastqc) and adapters and ambiguous bases were trimmed using Trimmomatic v0.39 (MINLEN: 75, SLIDINGWINDOW: 4:15, LEADING: 10, TRAILING: 10, AVGQUAL: 30). The trimmed RNA-seq reads were also assessed with FastQC before further analysis.

## 2. Files in this repository
  - **HMED_transdecoder.filtered.cds**
  - **HMED.mod.trinity.fasta**
