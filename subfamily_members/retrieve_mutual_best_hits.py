import os
import sys
import logging
import subprocess
import json
import tempfile

base_path = ""
json_files_path = base_path + "single_isoform_blastp_json_files/"
mutual_best_hits_files_path = base_path + "single_isoform_mutual_best_hits/"
fasta_files_path = base_path + "files_with_lineage_info/"


def retrieve_query_hits(query_tax_id,query_gene_id,subject_tax_id):

    query_json_file = json_files_path + query_tax_id + "_" + subject_tax_id + ".json"

    with open(query_json_file,'r') as outfile:


        query_json_dict = json.load(outfile)

    outfile.close()

    query_acc_no = list(query_json_dict[query_tax_id][query_gene_id].keys())[0]
    query_hit_tax_id = list(query_json_dict[query_tax_id][query_gene_id][query_acc_no].keys())[0]
    query_hit_gene_id = list(query_json_dict[query_tax_id][query_gene_id][query_acc_no][query_hit_tax_id].keys())[0]
    query_hit_acc_no = list(query_json_dict[query_tax_id][query_gene_id][query_acc_no][query_hit_tax_id][query_hit_gene_id].keys())[0]

    return query_acc_no, query_hit_acc_no, query_hit_gene_id, query_hit_tax_id



def retrieve_subject_hits(subject_tax_id,query_tax_id,query_hit_gene_id, query_hit_acc_no):

    subject_json_file = json_files_path + subject_tax_id + "_" + query_tax_id + ".json"
    
    with open(subject_json_file) as f:
        subject_json_dict = json.load(f)

    f.close()
    subject_hit_gene_id =  list(subject_json_dict[subject_tax_id][query_hit_gene_id][query_hit_acc_no][query_tax_id].keys())[0]
    subject_hit_acc_no = list(subject_json_dict[subject_tax_id][query_hit_gene_id][query_hit_acc_no][query_tax_id][subject_hit_gene_id].keys())[0]

    return subject_hit_acc_no



def get_fasta_dict(tax_id):

    fasta_dict = {}

    with open(fasta_files_path + tax_id + "_with_tax_info.fa",'r' ) as file:

        for line in file:
            if line.startswith(">"):
                header = line.strip()
                fasta_dict[header] = ''
            else:
                fasta_dict[header] += line.strip()
        
    file.close()

    return fasta_dict




def write_mutual_best_hits_to_fasta_file(query_tax_id,query_acc_no,subject_tax_id,query_hit_acc_no,query_gene_id):
    query_dict = get_fasta_dict(query_tax_id)
    subject_dict = get_fasta_dict(subject_tax_id)
    shorten_query_keys = [x.split("_geneid")[0] for x in query_dict.keys()]
    shorten_subject_keys = [y.split("_geneid")[0] for y in subject_dict.keys()]
    mutual_best_hits_file = mutual_best_hits_files_path + query_tax_id + "_" + query_gene_id + ".fa" 

    _ = open(mutual_best_hits_file, "a")
    _.close()

    with open(mutual_best_hits_file,'r') as f:
        file_content = f.read()
    f.close()

    with open(mutual_best_hits_file,'a') as f:
        if ">" + query_acc_no in shorten_query_keys:
            i = shorten_query_keys.index(">" + query_acc_no)
            if query_acc_no not in file_content:
                f.write(list(query_dict.keys())[i] + '\n' + query_dict[list(query_dict.keys())[i]] + '\n')

        if ">" + query_hit_acc_no in shorten_subject_keys:
            j = shorten_subject_keys.index(">" + query_hit_acc_no)
            if query_hit_acc_no not in file_content:
                f.write(list(subject_dict.keys())[j] + '\n' + subject_dict[list(subject_dict.keys())[j]]+'\n')
        
    f.close()
        
                

def retrieve_mutual_best_hits(query_tax_id,query_gene_id):
    fasta_files = os.listdir(fasta_files_path)
    for fasta_file in fasta_files:
        subject_tax_id = fasta_file.split("_")[0]

        if query_tax_id == subject_tax_id:
            continue
        
        query_acc_no, query_hit_acc_no, query_hit_gene_id ,query_hit_tax_id = retrieve_query_hits(query_tax_id,query_gene_id,subject_tax_id)
        if query_hit_gene_id == "NA" or query_hit_acc_no == "NA" or query_hit_tax_id == "NA":
            continue
        
        subject_hit_acc_no = retrieve_subject_hits(subject_tax_id,query_tax_id,query_hit_gene_id, query_hit_acc_no)
        if query_acc_no == subject_hit_acc_no:
            write_mutual_best_hits_to_fasta_file(query_tax_id,query_acc_no,subject_tax_id,query_hit_acc_no,query_gene_id)

       

if __name__ == "__main__":
    tax_id = sys.argv[1]
    gene_id = sys.argv[2]
    retrieve_mutual_best_hits(tax_id,gene_id)











        





