#!/usr/env/python

#This scripts gets the whole length of the sequence translated (taking frame into consideration), not only the fragment that has a hit with the database

from Bio import SeqIO

#Input files
blast_file="CMC-Chapaev-3-families_frame.blastout"
fasta_file="CMC-Chapaev-3-families.fa"

#Output file
pep_file="CMC-Chapaev-3-families_frame.pep"

seqs={}
#Read seqs
for seq in SeqIO.parse(fasta_file,'fasta'):
	#Add seqs to dictionary
	seqs[seq.id]=seq.seq


def extract_seq_and_translate(seq,frame):
	if frame<0:
		seq=seq.reverse_complement()
	if abs(frame)==2:
		seq=seq[1:]
	if abs(frame)==3:
		seq=seq[2:]
	return seq.translate(table='Standard', stop_symbol='*', to_stop=False)

peps={}
peps_length={}
#Read BLAST results
with open(blast_file,"rt") as blast_results:
	line=blast_results.readline().strip()
	while line:
		qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,frame=line.split("\t")
		if qseqid not in peps:
			pep_seq=extract_seq_and_translate(seqs[qseqid],int(frame))
			peps[qseqid]=pep_seq
			peps_length[qseqid]=length
		else:
			next
		line=blast_results.readline().strip()


with open(pep_file,"wt") as pep_out:
	for seq_id in peps.keys():
		if len(str(peps[seq_id]))>0 and int(peps_length[seq_id])>50:
			pep_out.write(">"+seq_id+"\n")
			pep_out.write(str(peps[seq_id])+"\n")
