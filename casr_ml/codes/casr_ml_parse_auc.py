import os
import sys
import csv


auc_dict = {}

for i in range(1001,1051):
    with open('/Users/aylin/Desktop/casr_artcile_ml_check_11_06_2023/with_pr_50_split_results_depth_2/casr_xgb_replication_'+str(i)+"_state.csv") as f:
        lines = f.readlines()
        for line in lines:
            print(line)
            if "train-auc" in line:
                train_auc = float(line.split("train-auc:")[1].split("\t")[0])
                test_auc = float(line.split("test-auc:")[1].split("\n")[0].split('"')[0])
                auc_dict[i] = {'train_auc':train_auc,'test_auc':test_auc}
    
print(auc_dict)

with open("/Users/aylin/Desktop/casr_artcile_ml_check_11_06_2023/with_pr_50_split_results_depth_2/casr_xgb_train_test_auc_2.csv","w") as f:
    auc_writer = csv.writer(f, delimiter=',')
    auc_writer.writerow(['replication', 'train_auc', 'test_auc'])
    for i in range(1001,1051):
        auc_writer.writerow([i,auc_dict[i]['train_auc'],auc_dict[i]['test_auc']])

f.close()

with open("/Users/aylin/Desktop/casr_artcile_ml_check_11_06_2023/with_pr_50_split_results_depth_2/casr_xgb_train_test_auc.csv","w") as f:
    auc_writer = csv.writer(f, delimiter=',')
    auc_writer.writerow(['replication', 'type', 'AUROC'])
    j = 1
    for i in range(1001,1051):
        auc_writer.writerow([j,'train',auc_dict[i]['train_auc']])
        auc_writer.writerow([j,'test',auc_dict[i]['test_auc']])
        j += 1

f.close()






pr_auc_dict = {}

for i in range(1001,1051):
    with open("/Users/aylin/Desktop/casr_artcile_ml_check_11_06_2023/with_pr_50_split_results_depth_2/casr_xgb_replication_" +str(i)+"_test_prauc.csv") as f:
        lines = f.readlines()
        for line in lines:
            #print(line)
            if '"1"' in line:
                test_pr_auc = float(line.split(",")[1].strip())
                pr_auc_dict[i] = {'test_pr_auc':test_pr_auc}
    with open("/Users/aylin/Desktop/casr_artcile_ml_check_11_06_2023/with_pr_50_split_results_depth_2/casr_xgb_replication_" +str(i)+"_train_prauc.csv") as f:
        lines = f.readlines()
        for line in lines:
            print(line)
            if '"1"' in line:
                train_pr_auc = float(line.split(",")[1].strip())
                pr_auc_dict[i].update({'train_pr_auc':train_pr_auc})
    
#print(pr_auc_dict)

with open("/Users/aylin/Desktop/casr_artcile_ml_check_11_06_2023/with_pr_50_split_results_depth_2/casr_xgb_train_test_prauc_2.csv","w") as f:
    auc_writer = csv.writer(f, delimiter=',')
    auc_writer.writerow(['replication', 'train_pr_auc', 'test_pr_auc'])
    for i in range(1001,1051):
        auc_writer.writerow([i,pr_auc_dict[i]['train_pr_auc'],pr_auc_dict[i]['test_pr_auc']])

f.close()

with open("/Users/aylin/Desktop/casr_artcile_ml_check_11_06_2023/with_pr_50_split_results_depth_2/casr_xgb_train_test_prauc.csv","w") as f:
    auc_writer = csv.writer(f, delimiter=',')
    auc_writer.writerow(['replication', 'type', 'AUPR'])
    j = 1
    for i in range(1001,1051):
        auc_writer.writerow([j,'train',pr_auc_dict[i]['train_pr_auc']])
        auc_writer.writerow([j,'test',pr_auc_dict[i]['test_pr_auc']])
        j += 1

f.close()


