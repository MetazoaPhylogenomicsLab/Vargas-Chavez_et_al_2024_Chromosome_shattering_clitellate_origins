#!/usr/bin/env/python

input_file="DEGs_components_num.txt"
composites_file="full_blastp_out_cleanNetwork_composites_Thu_Mar__7_07_42_58_2024/full_blastp_out_cleanNetwork.composites"
gene_hog_file="annelida_genenum_hog.txt"
gene_num_file="full_blastp.out.cleanNetwork.dico"

#Obtain gene names
gene_names={}
with open(gene_num_file, "rt") as convfile:
	line=convfile.readline().strip()
	while line:
		name,num=line.split()
		if not gene_names:
			gene_names={num:name}
		else:
			gene_names[num]=name
		line=convfile.readline().strip()

#Read list of components to analyze
components=[]
with open(input_file, "rt") as infile:
	line=infile.readline().strip()
	while line:
		if not components:
			components=[line]
		else:
			components.append(line)
		line=infile.readline().strip()

#Obtain which HOG each gene belongs to
gene_hog={}
hog_genes={}
with open(gene_hog_file, "rt") as convfile:
	line=convfile.readline().strip()
	while line:
		gene,hog=line.split(" ")
		if not gene_hog:
			gene_hog={"C"+gene:hog}
		else:
			gene_hog["C"+gene]=hog
		if not hog_genes:
			hog_genes={hog:[gene]}
		else:
			if hog not in hog_genes:
				hog_genes[hog]=[gene]
			else:
				hog_genes[hog].append(gene)
		line=convfile.readline().strip()

#Add genes from same HOG as components
components_new=[]
for component in components:
	genes=hog_genes[gene_hog["C"+component]]
	for gene in genes:
		if not components_new:
			components_new=[gene]
		else:
			components_new.append(gene)

components=components_new

#Obtain composite genes the input genes are components of
components_composites={}
with open(composites_file,"rt") as infile:
        line=infile.readline().strip()
        while line:
                if line[0]==">":
                        composite=line[1:]
                elif line[0]=="[":
                        line=infile.readline().strip()
                        continue
                else:
                        component=line.split("\t")[1]
                        if component not in components:
                                line=infile.readline().strip()
                                continue
                        else:
                                if not components_composites:
                                        components_composites={component:[composite]}
                                elif component not in components_composites:
                                        components_composites[component]=[composite]
                                else:
                                        components_composites[component].append(composite)
                line=infile.readline().strip()

#Print results
for component in components_composites.keys():
	for composite in components_composites[component]:
		print("Component "+str(gene_names[component])+" from "+str(gene_hog["C"+component])+"("+str(len(hog_genes[gene_hog["C"+component]]))+")"+": Composite "+str(gene_names[composite[1:]])+"("+str(gene_hog[composite])+"("+str(len(hog_genes[gene_hog[composite]]))+"))\n")
