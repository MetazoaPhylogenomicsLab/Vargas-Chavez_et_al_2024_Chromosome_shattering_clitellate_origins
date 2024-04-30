'''
2.2.count_table_from_pyham_output.py
script to obtain count table by HOG from OMA and pyHam outputs
wrote by Lisandra Benítez Álvarez on July 2023 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, Catalonia, Spain
'''

import os
from collections import Counter
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Obtain sequences count table from sequence composition files')
    parser.add_argument('-d', type=str, help='path to directory ')
    parser.add_argument('-o', type=str, help='output file')
    args = parser.parse_args()

    # Directory containing the txt files
    directory = args.d

    # Get a list of txt files in the directory
    txt_files = [filename for filename in os.listdir(directory) if filename.endswith('.txt')]

    # Initialize a dictionary to store count tables for each file
    count_tables = {}

    # Loop through each txt file
    for filename in txt_files:
        file_path = os.path.join(directory, filename)
        
        # Read the contents of the txt file
        with open(file_path, 'r') as file:
            elements = file.read().splitlines()
            elements = [item.split('.')[0] for item in elements]
        
        # Create a count table for the current file
        count_table = Counter(elements)
        count_tables[filename] = count_table

    # Get the unique elements across all files
    all_elements = set()
    for count_table in count_tables.values():
        all_elements.update(count_table.keys())
    unique_elements = sorted(all_elements)

    # Print the header row
    header_row = ['HOG'] + unique_elements

    # Print the count tables for each file
    for filename, count_table in count_tables.items():
        name = filename.replace("_seqs.txt", "")
        count_row = [name] + [str(count_table.get(element, 0)) for element in unique_elements]

    # Create and write the output to a TSV file
    output_file = args.o
    with open(output_file, 'w') as outfile:
        header_row = ['HOG'] + unique_elements
        outfile.write('\t'.join(header_row) + '\n')
        
        for filename, count_table in count_tables.items():
            name = filename.replace("_seqs.txt", "")
            count_row = [name] + [str(count_table.get(element, 0)) for element in unique_elements]
            outfile.write('\t'.join(count_row) + '\n')

    print(f"Count table saved to '{output_file}'.")

