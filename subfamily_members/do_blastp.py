import os
import sys
import logging
import subprocess
import json
import tempfile

def do_blastp(query, subject):
    base_path = ""
    blast_dbs_path = base_path + "blastdb_single_isoform/"
    fasta_files_for_blastp_path = base_path + "single_isoform/"
    json_file_path = base_path + "single_isoform_blastp_json_files/"
    jsons = []
    json_dict = {}
    fastaDict  = {}
    filein = open(fasta_files_for_blastp_path + query)
    json_files = os.listdir(json_file_path)
    logging.basicConfig(filename='json.log',level=logging.ERROR)
    for line in filein:
        if line.startswith('>'):
            header = line.strip()
            fastaDict[header] = ''
        else:
            fastaDict[header] += line.strip()

    for header in fastaDict.keys():
        #print(fastaDict)
        acc_no = header.split(">")[1].split("_geneid")[0]
        #print(acc_no)
        tax_id = header.split(">")[1].split("taxid")[1].split("_")[1]
        gene_id = header.split(">")[1].split("_geneid")[1].split("_")[1]

        fp = tempfile.NamedTemporaryFile(mode='r+')
        temp_out_file = tempfile.NamedTemporaryFile(mode='r+')
        
        fp.write(header +"\n" + fastaDict[header])
 
        fp.seek(0)
        #fp.read()
        
        bashCommand_tbl = "blastp -query " + fp.name+" -db " + blast_dbs_path + subject + " -outfmt 0 -out " + temp_out_file.name
        process_1 = subprocess.Popen(bashCommand_tbl.split(), stdout=subprocess.PIPE).wait()
        
        fp.close()
        
        best_saved = False
        with open(temp_out_file.name) as filein_:
            for line in filein_:
                if best_saved:
                    break

                if line.startswith(">"):
                    #print(line)
                    best_hit_acc_no = line.split(">")[1].split("_geneid")[0]
                    best_hit_tax_id = line.split("taxid_")[1].split("_")[0]
                    best_hit_gene_id = line.split("geneid_")[1].split("_")[0]
                    next(filein_)
                    next(filein_)
                    line = next(filein_)
                    score = line.split("Score = ")[1].split(",")[0].split("(")[1].split(")")[0]
                    e_value = line.split("Expect = ")[1].split(",")[0]
                    best_saved = True
                    #print(tax_id)
                    #print(gene_id)
                    #print(acc_no)
                

                    if query.split(".")[0].split("_")[0] +"_" +subject +".json" in json_files:
                        with open(json_file_path+query.split(".")[0].split("_")[0] +"_" +subject +".json", 'r') as f:
                            #print(f)

                            json_dict = json.load(f)
                            f.close()
                            #print(json_dict)
                


                    if tax_id in json_dict.keys():
                        if gene_id in json_dict[tax_id].keys():
                            if acc_no in json_dict[tax_id][gene_id].keys():
                                if best_hit_tax_id in json_dict[tax_id][gene_id][acc_no].keys():
                                    json_dict[tax_id][gene_id][acc_no][best_hit_tax_id].update(
                                        {
                                                best_hit_gene_id: {
                                                    best_hit_acc_no: {
                                                    "score":score,
                                                    "e_value":e_value
                                                    }
                                                }
                                            })
                                else:
                                    json_dict[tax_id][gene_id][acc_no].update(
                                        {
                                            best_hit_tax_id:{
                                                best_hit_gene_id: {
                                                    best_hit_acc_no: {
                                                    "score":score,
                                                    "e_value":e_value
                                                    }
                                                }
                                            }
                                        })
                            else:
                                json_dict[tax_id][gene_id].update({acc_no:{
                                                best_hit_tax_id: {
                                                    best_hit_gene_id: {
                                                        best_hit_acc_no: {
                                                        "score":score,
                                                        "e_value":e_value
                                                        }
                                                    }
                                                }
                                            }
                                })
                        else:
                            json_dict[tax_id].update({gene_id:{acc_no:{
                                                best_hit_tax_id: {
                                                    best_hit_gene_id: {
                                                        best_hit_acc_no: {
                                                        "score":score,
                                                        "e_value":e_value
                                                        }
                                                    }
                                                }
                                            }
                                }
                            })
            
                    
                    else:
                        
                        if tax_id in json_dict.keys():
                            if gene_id in json_dict[tax_id].keys():
                                if acc_no in json_dict[tax_id][gene_id].keys():
                                    if best_hit_tax_id in json_dict[tax_id][gene_id][acc_no].keys():
                                        json_dict[tax_id][gene_id][acc_no][best_hit_tax_id].update({
                                                    best_hit_gene_id: {
                                                        best_hit_acc_no: {
                                                        "score":score,
                                                        "e_value":e_value
                                                        }
                                                    }
                                            })
                                    else:
                                        json_dict[tax_id][gene_id][acc_no].update({best_hit_tax_id:{
                                                    best_hit_gene_id: {
                                                        best_hit_acc_no: {
                                                        "score":score,
                                                        "e_value":e_value
                                                        }
                                                    }
                                        }
                                            })
                                        
                                else:
                                    json_dict[tax_id][gene_id].update({acc_no:{
                                                    best_hit_tax_id: {
                                                        best_hit_gene_id: {
                                                            best_hit_acc_no: {
                                                            "score":score,
                                                            "e_value":e_value
                                                            }
                                                        }
                                                    }
                                                }
                                    })
                            else:
                                json_dict[tax_id].update({gene_id:{acc_no:{
                                                    best_hit_tax_id: {
                                                        best_hit_gene_id: {
                                                            best_hit_acc_no: {
                                                            "score":score,
                                                            "e_value":e_value
                                                            }
                                                        }
                                                    }
                                                }
                                    }
                                })

                        else:
                            json_dict = {
                                tax_id: {
                                    gene_id: {
                                        acc_no: {
                                            best_hit_tax_id: {
                                                best_hit_gene_id: {
                                                    best_hit_acc_no: {
                                                    "score":score,
                                                    "e_value":e_value
                                                    }
                                                }
                                            }
                                        }
                                    }    
                                }
                            }
                        
                    temp_out_file.close()
                
                    with open(json_file_path + query.split(".")[0].split("_")[0] +"_" +subject+ ".json", 'w') as outfile:
                        js = json.dumps(json_dict,sort_keys=True, indent=4, separators=(',', ': '))
                        outfile.write(js)
                    outfile.close()
                if line.startswith("***** No hits found *****"):
                    best_hit_acc_no = "NA"
                    best_hit_tax_id = "NA"
                    best_hit_gene_id = "NA"
                   
                    score = "NA"
                    e_value = "NA"
                    best_saved = True


                    if query.split(".")[0].split("_")[0] +"_" +subject +".json" in json_files:
                        with open(json_file_path+query.split(".")[0].split("_")[0] +"_" +subject +".json", 'r') as f:
                            #print(f)

                            json_dict = json.load(f)
                            f.close()
                            #print(json_dict)
                


                    if tax_id in json_dict.keys():
                        if gene_id in json_dict[tax_id].keys():
                            if acc_no in json_dict[tax_id][gene_id].keys():
                                if best_hit_tax_id in json_dict[tax_id][gene_id][acc_no].keys():
                                    json_dict[tax_id][gene_id][acc_no][best_hit_tax_id].update(
                                        {
                                                best_hit_gene_id: {
                                                    best_hit_acc_no: {
                                                    "score":score,
                                                    "e_value":e_value
                                                    }
                                                }
                                            })
                                else:
                                    json_dict[tax_id][gene_id][acc_no].update(
                                        {
                                            best_hit_tax_id:{
                                                best_hit_gene_id: {
                                                    best_hit_acc_no: {
                                                    "score":score,
                                                    "e_value":e_value
                                                    }
                                                }
                                            }
                                        })
                            else:
                                json_dict[tax_id][gene_id].update({acc_no:{
                                                best_hit_tax_id: {
                                                    best_hit_gene_id: {
                                                        best_hit_acc_no: {
                                                        "score":score,
                                                        "e_value":e_value
                                                        }
                                                    }
                                                }
                                            }
                                })
                        else:
                            json_dict[tax_id].update({gene_id:{acc_no:{
                                                best_hit_tax_id: {
                                                    best_hit_gene_id: {
                                                        best_hit_acc_no: {
                                                        "score":score,
                                                        "e_value":e_value
                                                        }
                                                    }
                                                }
                                            }
                                }
                            })
            
                    
                    else:
                        
                        if tax_id in json_dict.keys():
                            if gene_id in json_dict[tax_id].keys():
                                if acc_no in json_dict[tax_id][gene_id].keys():
                                    if best_hit_tax_id in json_dict[tax_id][gene_id][acc_no].keys():
                                        json_dict[tax_id][gene_id][acc_no][best_hit_tax_id].update({
                                                    best_hit_gene_id: {
                                                        best_hit_acc_no: {
                                                        "score":score,
                                                        "e_value":e_value
                                                        }
                                                    }
                                            })
                                    else:
                                        json_dict[tax_id][gene_id][acc_no].update({best_hit_tax_id:{
                                                    best_hit_gene_id: {
                                                        best_hit_acc_no: {
                                                        "score":score,
                                                        "e_value":e_value
                                                        }
                                                    }
                                        }
                                            })
                                        
                                else:
                                    json_dict[tax_id][gene_id].update({acc_no:{
                                                    best_hit_tax_id: {
                                                        best_hit_gene_id: {
                                                            best_hit_acc_no: {
                                                            "score":score,
                                                            "e_value":e_value
                                                            }
                                                        }
                                                    }
                                                }
                                    })
                            else:
                                json_dict[tax_id].update({gene_id:{acc_no:{
                                                    best_hit_tax_id: {
                                                        best_hit_gene_id: {
                                                            best_hit_acc_no: {
                                                            "score":score,
                                                            "e_value":e_value
                                                            }
                                                        }
                                                    }
                                                }
                                    }
                                })

                        else:
                            json_dict = {
                                tax_id: {
                                    gene_id: {
                                        acc_no: {
                                            best_hit_tax_id: {
                                                best_hit_gene_id: {
                                                    best_hit_acc_no: {
                                                    "score":score,
                                                    "e_value":e_value
                                                    }
                                                }
                                            }
                                        }
                                    }    
                                }
                            }
                        
                    temp_out_file.close()
                
                    with open(json_file_path + query.split(".")[0].split("_")[0] +"_" +subject+ ".json", 'w') as outfile:
                        js = json.dumps(json_dict,sort_keys=True, indent=4, separators=(',', ': '))
                        outfile.write(js)
                    outfile.close()
        
if __name__ == "__main__":
    q = sys.argv[1] 
    s = sys.argv[2]
    do_blastp(q,s)