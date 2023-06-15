import itertools
import os

GPCR_files_path = "GPCR_fasta_files/"
GPCR_files = os.listdir(GPCR_files_path)
single_isoforms_path = "single_isoform/" 




def parse_header(header):
    tax_id = header.split("taxid_")[1].split("_")[0]
    gene_id = header.split("geneid_")[1].split("_")[0]
    acc_no = header.split("geneid")[0]
    name = header.split("taxid_")[1].split("_")[1:]
    name = '_'.join(name)
    return {
        "gene_id":gene_id,
        "tax_id":tax_id,
        "acc_no":acc_no,
        "name":name
        }

def fasta_reader(filename):
    sequences = []
    seqDict = {}
    filein = open(GPCR_files_path + filename, 'r')
    # for line in filein:
    lines = filein.readlines()
    for i in range(0, len(lines)):
        match = False
        if lines[i].startswith('>'):
            seqDict = parse_header(lines[i])
            sequence = lines[i+1]
            seqDict["sequence"] = sequence
            #print(seqDict["name"])
            #print(len(seqDict["sequence"]))

            for seq in sequences:
                if seq["tax_id"] == seqDict["tax_id"] and seq["gene_id"] == seqDict["gene_id"]: 
                    if len(seqDict["sequence"]) > len(seq["sequence"]):
                        # print(len(seqDict["sequence"]))
                        # print(len(seq["sequence"]))
                        seq["acc_no"] = seqDict["acc_no"]
                        seq["name"] = seqDict["name"]
                        seq["sequence"] = seqDict["sequence"]
                    match = True
                    continue

            if not match:
                sequences.append(seqDict)
    #print(sequences)
    return sequences

for file in GPCR_files:
    
    if "_with_tax_info" not in file:
        continue 

    single_isoform_list = fasta_reader(file)

    single_isoform_fasta_file = open(single_isoforms_path + file ,"w")

    for i in single_isoform_list:
        single_isoform_fasta_file.write(i["acc_no"].strip() + "taxid_" + i["tax_id"].strip() +"_" + "geneid_" + i["gene_id"].strip()  + "_" + i["name"]+  i["sequence"].strip() + "\n")


    



