import os
import sys
import logging
import subprocess
import json
import operator
from fasta_dict import *
from internal_node_info import *
import ete3
from ete3 import Tree
from ete3 import NCBITaxa
ncbi = NCBITaxa()
base_path = ""
json_file_path = base_path + "hmmscan_json_file/"
subclass_file_path = base_path + "subclass_sequences_02/"
files_with_lineage_path = base_path + "single_isoform/"
tree_files_path = base_path + "subclass/"
trees_with_tax_path = base_path +"subclass_trees_with_most_specific_tax/"
def get_json_file(json_file):
    with open(json_file_path + json_file,'r') as outfile:
        json_dict = json.load(outfile)

    return json_dict

def get_sequences_with_highest_score(json_dict):
    
    subclass_dict = {}
    
    for taxid in json_dict.keys():
        subclass_dict[taxid] = {}
        for geneid in json_dict[taxid].keys():
            subclass_dict[taxid][geneid] = {}
            for acc_no in json_dict[taxid][geneid].keys():
                best_score = None
                best_e_value = None
                best_query = None
                best_lineage = None
                
                for query, value in json_dict[taxid][geneid][acc_no].items():
                    if best_score == None:
                        best_score = json_dict[taxid][geneid][acc_no][query]['score']
                    if best_e_value == None:
                        best_e_value = json_dict[taxid][geneid][acc_no][query]['e_value']

                    if best_query == None:
                        best_query = query

                    if best_lineage == None:
                        best_lineage = json_dict[taxid][geneid][acc_no][query]['lineage']
                        
                    if json_dict[taxid][geneid][acc_no][query]['e_value'] <= best_e_value and json_dict[taxid][geneid][acc_no][query]['score'] >= best_score:
                        best_score = json_dict[taxid][geneid][acc_no][query]['score']
                        best_e_value = json_dict[taxid][geneid][acc_no][query]['e_value']
                        best_query = query
                        best_lineage = json_dict[taxid][geneid][acc_no][query]['lineage']
                    
                        subclass_dict[taxid][geneid][acc_no] = {query: value}

                        

    return subclass_dict

def check_subclass(subclass_dict,tree_file,subclass,most_common_tax_rank):
    subclass_seq_dict = {}
    for taxid in subclass_dict.keys():
        for geneid in subclass_dict[taxid].keys():
            for acc_no in subclass_dict[taxid][geneid].keys():
                for query, value in subclass_dict[taxid][geneid][acc_no].items():
                    if query == subclass:
                        if int(most_common_tax_rank.split("_")[0]) in value['lineage']:
                            if taxid not in subclass_seq_dict.keys():
                                subclass_seq_dict[taxid]= []
                                subclass_seq_dict[taxid].append(acc_no)
                            else:
                                subclass_seq_dict[taxid].append(acc_no)

                                
    return subclass_seq_dict





def write_to_file(subclass_seq_dict,single_isoform_dict):
    for taxid,acc_nos in subclass_seq_dict.items():
        for acc_no in acc_nos:
            for header in single_isoform_dict.keys():
                if acc_no in header:
                    if "synaptonemal_complex" in header:
                        continue 
                    if "fer_1_like" in header:
                        continue
                    if "serine_protease_27_like" in header:
                        continue
                    if  "adhesion_G_protein"  in header:
                        continue
                    if "growth/differentiation" in header:
                        continue
                    with open(subclass_file_path + subclass + "_sequences.fa",'a') as f:
                        f.write(header + '\n' + single_isoform_dict[header]+ '\n')
                    f.close()


if __name__ == "__main__":
    json_dict =  get_json_file("sequences_against_GPCRdb.json")
    tree_file = sys.argv[1]
    subclass = sys.argv[2]
    single_isoform_dict = get_fasta_dict(files_with_lineage_path,"compiled_single_isoform_with_lineage_info.fa")
    most_common_tax_rank =  get_most_common_specific_tax_rank(tree_file)
    subclass_dict = get_sequences_with_highest_score(json_dict)
    subclass_seq_dict = check_subclass(subclass_dict,tree_file,subclass,most_common_tax_rank)
    write_to_file(subclass_seq_dict,single_isoform_dict)
