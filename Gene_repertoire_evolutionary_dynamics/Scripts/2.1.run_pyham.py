'''
2.1.run_pyham.py
script to analyse OMA's output with pyHam
wrote by Lisandra Benítez Álvarez on July 2023 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology (IBE-UPF), CSIC
'''
import pyham
import csv
import re
import xml.etree.ElementTree as ET
import sys
import os
import argparse
import os
import shutil

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Obtain ancestral genomes from OMA's HOGs using pyHam")
    parser.add_argument("-x", type=str, help="HierarchicalGroups.orthoxml", required=True)
    parser.add_argument("-t", type=str, help="species tree with node annotation in nwk format", required=True)
    parser.add_argument("-n", type=str, help="list of nodes to obtain the ancestral genomes. One node by line matching the annotation in the species tree", required=True)
    parser.add_argument("-o", type=str, help="output dir")
    args = parser.parse_args()

    #Function to save results. Save the results just as is shown in the screen with print
    def save_data(filename, result):                    
        with open(filename, "w+") as output_file:
            print(result, file=output_file)

    #Create output dir 
    out_dir = args.o
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    #Import HOGs in orthoxml format
    HOG_xml_path=args.x

    #Import tree in nwk format
    tree_path=args.t
    tree_str = pyham.utils.get_newick_string(tree_path, type="nwk")
    print("species tree loaded")
    print(tree_str)

    #Import node list
    node_list_path = args.n

    #Create pyHam object
    ham_analysis = pyham.Ham(tree_str, HOG_xml_path, use_internal_name=True)
    print("ham_analysis object created")

    # Obtain TREEPROFILE
    treeprofile = ham_analysis.create_tree_profile(outfile = out_dir + "/treeprofile.html")
    print("treeprofile saved")

    #Print TreeMap
    myTreeMap = ham_analysis.taxonomy.tree.copy(method="newick")
    print("TreeMap visualization")
    print(myTreeMap)

    # Print node names 
    print("Ancestral genomes name using newick names:")
    for ag in ham_analysis.taxonomy.internal_nodes:
        print("\t- {}".format(ag.name))

    #Get all HOGs
    top_hogs = ham_analysis.top_level_hogs
    #print(top_hogs)
    print("Number of HOGs in the input data")
    print(len(top_hogs))
    save_data(filename = out_dir + "/top_hogs_pyham.txt", result=top_hogs)

    top_hogs_list = list(top_hogs.keys())
    top_hogs_list = [key.replace(":", "") for key in top_hogs_list]
    with open(out_dir + "/top_hogs_list.txt" , 'w') as file:
        for hog in top_hogs_list:
            file.write(hog + '\n')
    
    # Get the ancestral genomes of interest
    with open(args.n, "r") as file:
        node_list = [line.strip() for line in file]

    for node in node_list:
        print(node)

        anc_node = ham_analysis.get_ancestral_genome_by_name(node)
        print("Number of genes at ancestral genome in " + node)
        print(anc_node.get_number_genes())
        anc_node_genes = anc_node.genes
        save_data(filename= out_dir + "/" + node + "_ancestral.txt", result=anc_node_genes)

    #Decomposing orthoxml into species composition lists
    #Parsing orthoxml dat
    prefix = "{http://orthoXML.org/2011/}"
    tree = ET.parse(HOG_xml_path)
    root = tree.getroot()
    ortho_data = {}

    #Save a file by hog with sequence composition
    dir = out_dir + "/HOG_seq_comp_lists"
    os.makedirs(dir)
    for ortho_group in root.findall(f".//{prefix}orthologGroup"):
        group_id = ortho_group.get("id")
        name = group_id.replace(":", "")
        with open(f"{dir}/{name}_seqs.txt", "w+") as output_file:
            for gene in ortho_group.findall(f".//{prefix}geneRef"):
                gene_id = gene.get("id")
                gene_id_to_find = gene_id
                gene = root.find(".//{http://orthoXML.org/2011/}gene[@id='" + gene_id_to_find + "']")
                prot_id = gene.get("protId")
                print(prot_id, file=output_file)

        
