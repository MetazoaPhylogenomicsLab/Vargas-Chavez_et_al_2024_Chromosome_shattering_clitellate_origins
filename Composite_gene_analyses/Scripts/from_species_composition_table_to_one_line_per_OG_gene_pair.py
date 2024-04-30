#!/usr/bin/env/python

'''
Script created by Gemma I. Martinez-Redondo to obtain from species composition table (OrthoFinder format) a file with a line per pair of OG-gene.
Usage:
        python from_species_composition_table_to_one_line_per_OG_gene_pair.py -i infile1.txt [-o outfile.txt]
'''

import argparse
from pathlib import Path

#Define parsing of the command
def parser():
	args = argparse.ArgumentParser(description='Obtain pairwise combinations given a list.')
	args.add_argument('-i', '--infile', required=True, help="Name of the species composition table input file.")
	args.add_argument('-o', '--outfile', required=False, help="Path to the output file. If not provided, 'input1.out' will be used as default.")
	args=args.parse_args()
	return args

#Obtain arguments
file = parser().infile
outfile = parser().outfile

#Use defaults if optional arguments not given
if not outfile:
	p=Path(infile1)
	extensions="".join(p.suffixes)
	file_name=str(p).replace(extensions, "")
	outfile=file_name+".out"

#Read file
og_genes={}
with open(file,"rt") as infile:
	#Skip first line
	line=infile.readline().strip()
	line=infile.readline().strip()
	while line:
		columns=line.split("\t")
		og=columns[0]
		genes=columns[1:]
		if og not in og_genes.keys():
			og_genes[og]=[gene.split(",") for gene in genes]
			og_genes[og]=[k for i in og_genes[og] for k in i]
		else:
			og_genes[og].append(genes)
		line=infile.readline().strip()



#Print output file
with open(outfile,"wt") as out:
	for og in og_genes.keys():
		for gene in og_genes[og]:
			if gene != "-":
				out.write(gene.replace(" ","")+"\t"+og+"\n")
