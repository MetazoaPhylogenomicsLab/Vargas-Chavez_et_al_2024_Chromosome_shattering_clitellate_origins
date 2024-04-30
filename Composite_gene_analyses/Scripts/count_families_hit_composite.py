#!/usr/bin/env/python
composition_table={}

with open("HOG00040255_63_components.txt", "rt") as infile:
	line=infile.readline().strip()
	composite=line[1:]
	composition_table={composite:{}}
	line=infile.readline().strip()
	while line:
		if line[0]==">":
			composite=line[1:]
			composition_table[composite]={}
		elif line[0]=="[":
			line=infile.readline().strip()
			continue
		else:
			family=line.split(" ")[0]
			if family not in composition_table[composite]:
				composition_table[composite][family]=1
			else:
				composition_table[composite][family]+=1
		line=infile.readline().strip()

for composite in composition_table:
	print(str(composite)+": "+str(composition_table[composite]))
