#!/usr/bin/python

'''
Script created by Gemma I. Martinez-Redondo to obtain pairwise combinations given a list of items.
Usage:
        python obtain_pairwise_combinations_from_list.py -i infile.txt -o outfile.txt
'''

#Import required modules
import argparse, itertools
from pathlib import Path

#Define parsing of the command
def parser():
    args = argparse.ArgumentParser(description='Obtain pairwise combinations given a list.')
    args.add_argument('-i', '--infile', required=True, help="Name of the input file. It must have one line per item.")
    args.add_argument('-o', '--outfile', required=False, help="Name of the output file. If not provided, 'inputfilename.out' will be used as default.")
    args=args.parse_args()
    return args

#Obtain output file name
infile = parser().infile
outfile = parser().outfile

if not outfile:
    p=Path(infile)
    extensions="".join(p.suffixes)
    file_name=str(p).replace(extensions, "")
    outfile=file_name+".out"

items=set()
with open(infile,"rt") as in_file:
	line=in_file.readline().strip()
	while line:
		items.add(line)
		line=in_file.readline().strip()

combinations=itertools.product(items,repeat=2)

with open(outfile,"wt") as out_file:
	for combination in combinations:
		out_file.write(combination[0]+"\t"+combination[1]+"\n")
