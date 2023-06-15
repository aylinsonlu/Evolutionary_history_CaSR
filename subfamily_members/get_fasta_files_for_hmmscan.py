import os
hmmsearch_results_path= "/media/disk02/GPCR/hmmsearch_results/"
proteome_path= "/media/disk02/GPCR/proteome/"
fasta_files_for_hmmscan_path="/media/disk02/GPCR/fasta_files_for_hmmscan/"

files_in_hmmsearch= os.listdir(hmmsearch_results_path)
files_in_proteome= os.listdir(proteome_path)

acc_numbers=[]

for file in files_in_hmmsearch:
    print(file)
    if file + ".fa" in fasta_files_for_hmmscan_path:
        continue
    f=open (hmmsearch_results_path + file)
    for line in f:
        words=line.split()
        if "threshold" in words:
            break
        for word in words:
            if word.find("ref")>=0:
                if word not in acc_numbers:
                    acc_numbers.append(word)

    f= open (proteome_path + file.split(".out",1)[0])
    content=f.read()
    sequences=content.split(">")
    for acc_number in acc_numbers:
        for sequence in sequences:
            if acc_number in sequence:
                new_file=open(fasta_files_for_hmmscan_path + file.split(".out",1)[0] + ".fa","a")
                new_file.write(">" + sequence)
                print (file.split(".out",1)[0] + ".fa")
