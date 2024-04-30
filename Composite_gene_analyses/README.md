# Composite gene analyses

## 1. Description
To identify HOGs that arose through fusion or fission the following strategy was followed. An all vs all BLASTp of the proteomes was performed to create a sequence similarity network with the cleanBlastp command in CompositeSearch. The correspondence between the filtered HOGs and the genes was used, together with the output of cleanBlastp, as input for CompositeSearch to identify composite and component genes. Parameter values suggested in the tutorial were used (E-value of 1e-5, 30% identity, 80% coverage, maximum overlap of 20). A HOG was considered to be a composite HOG when more than half of the genes belonging to that HOG were composites. A composite HOG was considered to have originated from a fusion when most of its component HOGs were inferred to have originated before the origin of the composite HOG. Similarly, a fission origin is inferred when the age of origin of the component HOGs is younger than the age of origin of the component HOG. In the case of component HOGs, it was considered part of a fusion when all the genes in the HOG were components, and most of the composite genes these genes were components of had an origin younger than the origin of the component HOG. Contrarily, when the origin of the composite HOG predated the origin of the component HOG, it was considered to have originated from fission.
      
## 2. Commands:
Obtain list of annelida species pairwise combinations:
```
python obtain_pairwise_combinations_from_list.py -i annelida_species.txt -o annelida_pairs.txt
```
Execute all-vs-all BLASTp:
```
FILE=$WD/annelida_pairs.txt

SP1=$(awk -v num_line=$SLURM_ARRAY_TASK_ID '{if(NR==num_line) print $1}' $FILE)
SP2=$(awk -v num_line=$SLURM_ARRAY_TASK_ID '{if(NR==num_line) print $2}' $FILE)

IN_PATH=Annelida_project/final_annelida_pep_files

if [ -f $WD/blastp_results/${SP1}_${SP2}_blastp.out ]; then
	echo "$WD/blastp_results/${SP1}_${SP2}_blastp.out exists already"
else
	blastp -subject $IN_PATH/${SP1}.fa -query $IN_PATH/${SP2}.fa -seg yes -soft_masking true -max_target_seqs 5000 -outfmt "6 qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen" -out $WD/blastp_results/${SP1}_${SP2}_blastp.out
fi
```
Check that BLAST has executed correctly:
```
FILE=$WD/annelida_pairs.txt

while read -r LINE
do
	SP1=$(echo $LINE | awk '{print $1}')
	SP2=$(echo $LINE | awk '{print $2}')
	if [[ -f $WD/blastp_results/${SP1}_${SP2}_blastp.out ]]; then
		echo "${SP1}_${SP2}_blastp.out" >> blastp_done.txt
	else
		echo "${SP1}_${SP2}_blastp.out" >> blastp_missing.txt
	fi
done < $FILE
```
Combine BLASTp results:
```
cat $WD/blastp_results/* > $WD/full_blastp.out
```
Execute CompositeSearch:
```
#CleanBlastp
cleanblastp -i $WD/full_blastp.out -n 1

#Obtain gene-gene family file needed for CompositeSearch
python $WD/scripts/from_species_composition_table_to_one_line_per_OG_gene_pair.py -i $WD/species_composition_table_FILTERED_HOGs.tsv -o $WD/annelida_gene_hog.txt
awk 'FNR==NR { id[$1]=$2; next } { split($1,a,"\t"); if (a[1] in id) $1=id[a[1]]; print }' $WD/full_blastp.out.cleanNetwork.dico $WD/annelida_gene_hog.txt > $WD/annelida_genenum_hog.txt

#CompositeSearch (gene-gene family format: gene num in output of cleanblastp + "\t" + gene family
compositeSearch -i $WD/full_blastp.out.cleanNetwork -n $WD/full_blastp.out.cleanNetwork.genes -m composites -e 1e-05 -p 30 -c 80 -l 20 -t 24 -f $WD/annelida_genenum_hog.txt
```
Obtain data for Supplementary table:
```
#Get HOG names for both composites
grep -f DEGs_composites.txt species_composition_table_FILTERED_HOGs.tsv

#Get genes from same HOGs to see if they are composites
#EAND.evm.TU.Chr10.1856 (237136)
grep HOG00040281_63 annelida_genenum_hog.txt | cut -f1 -d" " > HOG00040281_63_genes.txt
grep -f HOG00040281_63_genes.txt full_blastp_out_cleanNetwork_composites_Thu_Mar__7_07_42_58_2024/full_blastp_out_cleanNetwork.compositesinfo
#Result: only C237136
#EAND.evm.TU.Chr03.3017 (221222)
grep HOG00040255_63 annelida_genenum_hog.txt | cut -f1 -d" " > HOG00040255_63_genes.txt
grep -f HOG00040255_63_genes.txt full_blastp_out_cleanNetwork_composites_Thu_Mar__7_07_42_58_2024/full_blastp_out_cleanNetwork.compositesinfo
# Result: C243566,C873971,C247941,C377041,C94033,C102461,C799630,C559357,C616388,C886794,C221222,C383931,C575467 (13/21 are composites)
# Copy in file with HOG00040255_63_composites.txt and without C HOG00040255_63_composites_num.txt
# Check which genes they are
grep -wf HOG00040255_63_composites_num.txt full_blastp.out.cleanNetwork.dico
# 2/2 MVUL, 1/1 EAND, 2/2 CMAT, 1/2 in NNAJ, 0/1 WPIG, 1/2 HNIP, 1/1 HMED, 0/2 AUJA, 1/1 PHRE, 1/1 PVOL, 2/2 ECRY, 0/1 HMAN, 0/2 TSTR, 1/1 RAND
#Check if these are components:
grep -wf HOG00040255_63_composites_num.txt full_blastp_out_cleanNetwork_composites_Thu_Mar__7_07_42_58_2024/full_blastp_out_cleanNetwork.composites
#These composite genes are components of C224234
grep -w 224234 full_blastp.out.cleanNetwork.dico #components of another EAND gene EAND.evm.TU.Chr04.2898 (HOG00103199_73) #Fusion in EAND

#Know the HOG corresponding to the family (numbers without F)
grep -w 24113 full_blastp_out_cleanNetwork_composites_Thu_Mar__7_07_42_58_2024/family.info # HOG00063769_68
grep -w 48213 full_blastp_out_cleanNetwork_composites_Thu_Mar__7_07_42_58_2024/family.info # HOG00071544_68
grep -w 6107 full_blastp_out_cleanNetwork_composites_Thu_Mar__7_07_42_58_2024/family.info # HOG00038038_62
grep -w 35306 full_blastp_out_cleanNetwork_composites_Thu_Mar__7_07_42_58_2024/family.info # HOG00010934_35
grep -w 50429 full_blastp_out_cleanNetwork_composites_Thu_Mar__7_07_42_58_2024/family.info # HOG00038195_62
grep -w 44347 full_blastp_out_cleanNetwork_composites_Thu_Mar__7_07_42_58_2024/family.info # HOG00002981_16

#Get number of composite HOGs each DEG is component of.
grep HOG Components_EAND_and_others_HOG_info_results.txt | cut -f3 -d"(" | sort | uniq -c | sort -rn
```

## 3. Files in this repository
  - **Files**:
    - **full_blastp.out.tar.gz**: CompositeSearch input
	- **full_blastp.out.cleanNetwork.tar.gz, full_blastp.out.cleanNetwork.dico.tar.gz** and **full_blastp.out.cleanNetwork.genes.tar.gz**: CompositeSearch preprocessing outputs
	- **full_blastp_out_cleanNetwork_composites_Thu_Mar__7_07_42_58_2024.tar.gz**: folder with outputs of CompositeSearch
    - **species_composition_table_FILTERED_HOGs.tsv**: Species composition table from filtered OMA results
    - **DEG_EAND_compositeSearch_num.txt**:  DEGs number assigned by CompositeSearch.
    - **DEG_EAND_conversion.txt**: Conversion of DEGs in EAND with numbers assigned by CompositeSearch.
    - **DEG_EAND.txt**: List of DEGs in EAND.
    - **DEGs_components_num.txt**: DEGs (number assigned by CompositeSearch) that are components.
    - **DEGs_components.txt**: DEGs that are components.
    - **DEGs_composites.txt**: DEGs that are composites.
    - **HOG00040255_63_components.txt**: components from HOG00040255.
    - **HOG00040255_63_composites_num.txt**: composite genes from HOG00040255 (number assigned by CompositeSearch).
    - **HOG00040255_63_composites.txt**: composite genes from HOG00040255 (number assigned by CompositeSearch preceded by a C.
    - **HOG00040255_63_genes.txt**: genes from HOG00040255 (number assigned by CompositeSearch).
    - **HOG00040281_63_genes.txt**: genes from HOG00040281 (number assigned by CompositeSearch).
    - **annelida_genenum_hog.txt**: gene-HOG file.
    - **Components_EAND_and_others_HOG_info_results**: information on composites genes are components of (including HOG those composite/component genes belong to and number of genes in each HOG).
    - **annelida_pairs.txt**: annelida species pairwise combinations.
  - **Scripts**:
    - **from_species_composition_table_to_one_line_per_OG_gene_pair.py**: script to obtain from species composition table (OrthoFinder format) a file with a line per pair of OG-gene.
    - **obtain_pairwise_combinations_from_list.py**: script to obtain pairwise combinations given a list.
    - **obtain_composites_given_list_of_components.py**: script to obtain composites given list of components. Output of this script was saved in **Components_EAND_and_others_HOG_info_results**.
    - **count_families_hit_composite.py**: script to count the number of hits to each family for each composite.
  
