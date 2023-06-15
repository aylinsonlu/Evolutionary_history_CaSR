import os
import sys
import operator
hmm_scan_results_path = "/media/disk02/GPCR/hmm_scan_results/"
hmm_scan_results= os.listdir(hmm_scan_results_path)
filename = "/media/disk02/GPCR/compiled_fasta_files/compiled_fasta_file_with_geneid_taxid.fa"
fastaDict = {}
acc_no_to_result_dict = {}
def getaccno(file):
    for line in f:
        words = line.split()
        if words[0] == "7tm_3":
            acc_no = words[3].split("ref|",1)[1][:-1]
            score = float(words[13])
            align_start = int(words[17])
            align_stop = int(words[18])

            if acc_no not in acc_no_to_result_dict.keys():
                acc_no_to_result_dict[acc_no] = {"score":score, "align_start":align_start, "align_stop":align_stop }
            else:
                acc_no_to_result_dict[acc_no]["score"] = max([acc_no_to_result_dict[acc_no]["score"], score])
                acc_no_to_result_dict[acc_no]["align_start"] = align_start
                acc_no_to_result_dict[acc_no]["align_stop"] = align_stop

    return acc_no_to_result_dict

for file in hmm_scan_results:
    if not file.endswith(".out"):
        continue
    f=open (hmm_scan_results_path + file )
    getaccno(f)
print(acc_no_to_result_dict)

filein = open(filename, 'r')
for line in filein:
    if line.startswith('>'):
        header = line[1:].strip()
        fastaDict[header] = ''
    else:
        fastaDict[header] += line.strip()
        
compiled_fasta_file_7tm_3_domains = open("/media/disk02/GPCR/compiled_fasta_files/compiled_fasta_file_7tm_3_domains.fa", 'w')
for header in fastaDict.keys():
    for acc_no,datadict in acc_no_to_result_dict.items():
        if acc_no in header:
            new_sequence = fastaDict[header][acc_no_to_result_dict[acc_no]["align_start"]-1:acc_no_to_result_dict[acc_no]["align_stop"]]
            compiled_fasta_file_7tm_3_domains.write('>' + header + '\n' + new_sequence +'\n')



    