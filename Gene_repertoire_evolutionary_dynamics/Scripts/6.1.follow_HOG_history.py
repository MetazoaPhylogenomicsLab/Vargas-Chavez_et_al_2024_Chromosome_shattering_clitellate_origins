'''
6.1.follow_HOG_history.py
script to obtain the ham graph and the tree profile for a list of HOGs
wrote by Lisandra Benítez Álvarez on August 2023 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, Catalonia, Spain
'''
import pyham
import sys
import os
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='obtain the ham graph and the tree profile for a list of HOGs')
    parser.add_argument('-d', type=str, help='path to directory containing HOGs lists. The files must have the extension .HOGs_list')
    parser.add_argument('-o', type=str, help='output dir')
    parser.add_argument('-x', type=str, help='path to orthoxml file')
    parser.add_argument('-t', type=str, help='path to tree file')
    args = parser.parse_args()

    def save_data(filename, result):                    
        with open(filename, "w+") as output_file:
            print(result, file=output_file)

    # Directory containing the list of files
    directory = args.d

    #Creating output directory
    out_dir = args.o
    os.makedirs(out_dir)

    #Import HOGs in orthoxml format
    HOG_xml_path=args.x

    #Import tree in nwk format
    tree_path=args.t
    tree_str = pyham.utils.get_newick_string(tree_path, type="nwk")
    print("species tree loaded")
    print(tree_str)

    #Create pyHam object
    ham_analysis = pyham.Ham(tree_str, HOG_xml_path, use_internal_name=True)
    print("ham_analysis object created")

    # Get a list of all text files in the directory
    file_names = [file for file in os.listdir(directory) if file.endswith('.HOGs_list')]
    print(file_names)

    # Process each list file
    for file_path in file_names:
        exp = file_path.replace(".HOGs_list", "")
        exp_dir = out_dir + "/" + exp
        os.makedirs(exp_dir)
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                line = line.replace('\n', '')
                HOG = line.replace("HOG", "HOG:")
                hog_ham =  ham_analysis.get_hog_by_id(HOG)
                hog_ham_desc = hog_ham.get_all_descendant_hogs()
                save_data(filename=exp_dir + '/' + line + '_descendents.txt', result=hog_ham_desc)
                ham_analysis.create_iHam(hog=hog_ham, outfile=exp_dir + '/' + line + '_ham.html'.format(hog_ham.hog_id))
                ham_analysis.create_tree_profile(hog=hog_ham, outfile=exp_dir + '/' + line + '_tp.html')
    