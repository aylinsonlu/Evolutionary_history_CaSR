import os
import sys
import operator

base_path = "/cta/users/aylinbircan/"
tm_domain_compiled_file_path = base_path + "GPCR/compiled_fasta_files/"
tm_domain_compiled_file = "compiled_fasta_file_7tm_3_domains.fa"
mutual_hits_path = base_path + "single_isoform_mutual_best_hits_with_9606/"
mutual_hits = os.listdir(mutual_hits_path)
mutual_hit_with_9606_7tm3_domain_path = base_path + "single_isoform_mutual_best_hits_with_9606_7tm3_domain/"

def get_fasta_dict(file_path,file):

    fasta_dict = {}

    with open(file_path + file,'r' ) as file:

        for line in file:
            if line.startswith(">"):
                header = line.strip()
                fasta_dict[header] = ''
            else:
                fasta_dict[header] += line.strip()
        
    file.close()

    return fasta_dict

def take_7tm3_domain(mutual_hit):

    tm3_compiled_dict = get_fasta_dict(tm_domain_compiled_file_path,tm_domain_compiled_file)
    shorten_compiled_keys = [y.split("ref|")[1].split("|")[0].replace(".","_") for y in tm3_compiled_dict.keys()] 
    mutual_hit_dict = get_fasta_dict(mutual_hits_path, mutual_hit)
    shorten_mutual_keys = [x.split("_geneid")[0].split(">")[1] for x in mutual_hit_dict.keys()]

    for acc_no in shorten_mutual_keys:
        i = shorten_mutual_keys.index(acc_no)

        if acc_no in shorten_compiled_keys:
            j = shorten_compiled_keys.index(acc_no)

            with open(mutual_hit_with_9606_7tm3_domain_path + mutual_hit,'a') as f:
                f.write(list(mutual_hit_dict.keys())[i] + '\n' + tm3_compiled_dict[list(tm3_compiled_dict.keys())[j]] + '\n')

            f.close()
        else:
            print(acc_no)

for mutual_hit in mutual_hits:
    take_7tm3_domain(mutual_hit)
