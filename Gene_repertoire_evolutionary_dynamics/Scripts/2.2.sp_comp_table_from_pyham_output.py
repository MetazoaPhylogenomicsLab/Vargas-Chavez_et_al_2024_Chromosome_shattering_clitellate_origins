'''
2.2.sp_comp_table_from_pyham_output.py
script to obtain species composition table by HOG from OMA and pyHam outputs
wrote by Lisandra Benítez Álvarez on July 2023 at Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, Catalonia, Spain
'''

import os
import pandas as pd
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Obtain sequences composition table from sequence composition files')
    parser.add_argument('-d', type=str, help='path to directory ')
    parser.add_argument('-o', type=str, help='output file')
    args = parser.parse_args()

    # Function to process a single list file
    def process_list_file(file_path):
        with open(file_path, 'r') as file:
            lines = file.readlines()
        
        names = {}
        for line in lines:
            #parts = line.strip().split('_') las seqs de Carlos tienen otro formato en el nombre
            parts = line.strip().split('.')
            if len(parts) > 1:
                name = parts[0]
                value = line.strip()
                if name in names:
                    names[name].append(value)
                else:
                    names[name] = [value]
        
        return names

    # Directory containing the list files
    directory = args.d

    # Get a list of all text files in the directory
    file_names = [file for file in os.listdir(directory) if file.endswith('.txt')]

    # Process each list file
    table = {}
    for file_name in file_names:
        file_path = os.path.join(directory, file_name)
        names = process_list_file(file_path)
        table[file_name] = names

    df = pd.DataFrame(table).T
    df.index.name = 'name'
    df.fillna('-', inplace=True)
    df.index = df.index.str.replace('_seqs.txt', '', regex=True)

    def process_cell(cell):
        if cell == '-':
            return '-'
        return ', '.join(item.strip('[]') for item in cell)

    df2 = df.applymap(process_cell)

    # Write the table to a TSV file
    output_file = args.o
    df2.to_csv(output_file, sep='\t', index=True)

    print(f"Output table saved to {output_file}")
