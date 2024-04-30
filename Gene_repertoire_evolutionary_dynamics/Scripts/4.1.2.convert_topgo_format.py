#!/usr/bin/python
'''
Script created by Gemma I. Martinez-Redondo to convert GOPredSim output format to topGO input format
Usage:
        python convert_topgo_format.py -a path_to_gopredsim_annotations -p prefix -o output_file.txt
'''

#Import required modules
import argparse
import os

#Define parsing of the command
def parser():
    args = argparse.ArgumentParser(description='Convert GOPredSim output format to topGO input format.')
    args.add_argument('-a', '--go_annotations', required=True, help="GOPredSim output annotations path. Take into account that the file names will be gopredsim_prefix_prott5_1_[bpo|mfo|cco].txt.")
    args.add_argument('-o', '--output_file', required=True, help="Output file name for new annotation.")
    args.add_argument('-p', '--prefix', required=True, help="Name to use as prefix to add to GOPredSim output folders.")
    args=args.parse_args()
    return args

#Obtain arguments and check their content
args=parser()

if not args.go_annotations:
    print("Path to output GOPredSim GO annotations must be provided")

if not args.output_file:
    print("Output file name must be provided")

if not args.prefix:
    print("Prefix must be provided")

#Read GO annotations and convert to topGO format
result = dict()
for category in ["bpo", "cco", "mfo"]:
    f = open(args.go_annotations + '/' + 'gopredsim_' + args.prefix + '_prott5_1_' + category + '.txt', 'r')
    for line in f:
        gene = line.split()[0]
        if gene not in result:
            result[gene] = set()
            result[gene].add(line.split()[1])
        else:
            result[gene].add(line.split()[1])

with open(args.output_file,"wt") as f_out:
    for k, v in result.items():
        f_out.write(f"{k}\t{', '.join(v)}\n")
