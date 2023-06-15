import sys
import os
import csv
import pandas as pd
import json
from fasta_dict import *
import operator
from Bio import AlignIO
from Bio.Align import AlignInfo
import pandas as pd


def get_aa_dict(msa_dict,position):
    aa_dict = {}
    for k,v in msa_dict.items():
        aa = v[position]
        if not aa in aa_dict.keys():
            aa_dict[aa] = 1
        else:
            aa_dict[aa] += 1
    if not "-" in aa_dict.keys():
        aa_dict["-"] = 0
    return aa_dict

def compute_identity_score(msa_dict,gap_treshold):
    seq_number = len(list(msa_dict.keys()))
    aln_lenght = len(list(msa_dict.values())[0])
    identity_dict = {}
    for i in range(aln_lenght):
        aa_dict = get_aa_dict(msa_dict,i)
        gap_prop = aa_dict["-"]/seq_number
        max_aa = max(aa_dict.items(), key=operator.itemgetter(1))[0]
        if max_aa == "-":
            identity_score = 0
            identity_dict[i+1] = {'identity_score':identity_score,'consensus_aa':max_aa}
        elif gap_prop <= gap_treshold:
            identity_score = aa_dict[max_aa]/(seq_number-aa_dict["-"])
            identity_dict[i+1] = {'identity_score':identity_score,'consensus_aa':max_aa}
        elif gap_prop > gap_treshold:
            identity_score = aa_dict[max_aa]/seq_number
            identity_dict[i+1] = {'identity_score':identity_score,'consensus_aa':max_aa}
    #print(identity_dict[725])
    return identity_dict

def write_csv(score_file,casr_fasta_file,casr_like_fasta_file,gprc6a_fasta_file,tas1r1_fasta_file,tas1r2_fasta_file,tas1r3_fasta_file):
    casr_dict = get_fasta_dict(casr_fasta_file)
    casr_identity = compute_identity_score(casr_dict,0.3)
    casr_like_dict = get_fasta_dict(casr_like_fasta_file)
    casr_like_identity = compute_identity_score(casr_like_dict,0.3)
    gprc6a_dict = get_fasta_dict(gprc6a_fasta_file)
    gprc6a_identity = compute_identity_score(gprc6a_dict,0.3)
    tas1r1_dict = get_fasta_dict(tas1r1_fasta_file)
    tas1r1_identity = compute_identity_score(tas1r1_dict,0.3)
    tas1r2_dict = get_fasta_dict(tas1r2_fasta_file)
    tas1r2_identity = compute_identity_score(tas1r2_dict,0.3)
    tas1r3_dict = get_fasta_dict(tas1r3_fasta_file)
    tas1r3_identity = compute_identity_score(tas1r3_dict,0.3)

    with open(score_file,'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow (["position","CaSR","CaSR_aa",
        "CaSR_like","CaSR_like_aa",           
        "GPRC6A","GPRC6A_aa",
        "TAS1R1","TAS1R1_aa",
        "TAS1R2","TAS1R2_aa",
        "TAS1R3","TAS1R3_aa"])
        for i in range(1,len(list(casr_identity.keys()))+1):
            print(i)
            writer.writerow([i,
            casr_identity[i]['identity_score'],casr_identity[i]['consensus_aa'],
            casr_like_identity[i]['identity_score'],casr_like_identity[i]['consensus_aa'],
            gprc6a_identity[i]['identity_score'],gprc6a_identity[i]['consensus_aa'],
            tas1r1_identity[i]['identity_score'],tas1r1_identity[i]['consensus_aa'],
            tas1r2_identity[i]['identity_score'],tas1r2_identity[i]['consensus_aa'],
            tas1r3_identity[i]['identity_score'],tas1r3_identity[i]['consensus_aa']])
    f.close()


write_csv("CASR.csv")


    