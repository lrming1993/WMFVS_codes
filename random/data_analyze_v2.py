import os
import numpy as np
import random

database_dict = {}

'''
CCGD_genes = []
census_genes = []
intogen_genes = []
oncogene_genes = []
TSG_genes = []
NCG_genes = []
CGN_genes = []
'''

# random_MFVS = []
gene_dict = {}
file_name_list = []
gene_lists = []
random_common_genes = []
random_half_common_genes = []
path = "random_result"

'''
gurobi_MFVS_genes = []
gurobi_WMFVS_genes = []
gurobi_WFVS_genes = []
'''

gene_weight = []
top_463_genes = []


for dirpath, dirnames, filenames in os.walk(path):
    for i in filenames:
        file_name_list.append(path + "/" + i)

for file_name in file_name_list:
    with open(file_name) as f0:
        temp = []
        for n, i in enumerate(f0):
            if n == 0:
                continue
            i = i.split()
            temp.append(i[0])
            if i[0] in gene_dict:
                gene_dict[i[0]] += 1
            else:
                gene_dict[i[0]] = 1
        gene_lists.append(temp)

for i in gene_dict:
    if gene_dict[i] == len(gene_lists):
        random_common_genes.append(i)
random_half_common_genes = list(i for i in gene_dict if gene_dict[i] >= 0.495 * len(gene_lists))

database_dict["random_common"] = random_common_genes
database_dict["random_half_common"] = random_half_common_genes

print("Got {} random files.".format(len(file_name_list)))
print("Got {} common genes.".format(len(database_dict["random_common"])))
print("Got {} genes appearing in 49.5% MFVSs".format(len(database_dict["random_half_common"])))

with open("gene_weight.txt") as f1:
    for n, i in enumerate(f1):
        i = i.split()
        gene_weight.append((i[0], float(i[-1])))

gene_weight = sorted(gene_weight, key=lambda x: x[-1])
top_463_genes = list(i[0] for i in gene_weight[-463:])
database_dict["top_463"] = top_463_genes
all_gene_list = list(i[0] for i in gene_weight)
database_dict["all"] = all_gene_list

print("Size of all genes:", len(all_gene_list))

def open_gene_file(database_name):
    database_dict[database_name] = []
    with open("database/{}_genes.txt".format(database_name)) as f1:
        for n, i in enumerate(f1):
            i = i.split()
            database_dict[database_name].append(i[0])
    print("Size of {}_genes:".format(database_name), len(database_dict[database_name]))


open_gene_file("CCGD")
open_gene_file("census")
open_gene_file("intogen")
open_gene_file("oncogene")
open_gene_file("TSG")
open_gene_file("NCG")

print("---------------")


def analyse_both(d1, d2):
    temp = list(i for i in database_dict[d1] if i in database_dict[d2])
    print("Both in {}_genes and {}_genes: {}".format(d1, d2, len(temp)))

analyse_both("all", "CCGD")
analyse_both("all", "TSG")
analyse_both("all", "oncogene")
analyse_both("all", "intogen")
analyse_both("all", "census")
analyse_both("all", "NCG")

print("---------------")

analyse_both("random_common", "CCGD")
analyse_both("random_common", "TSG")
analyse_both("random_common", "oncogene")
analyse_both("random_common", "intogen")
analyse_both("random_common", "census")
analyse_both("random_common", "NCG")

print("---------------")

analyse_both("top_463", "CCGD")
analyse_both("top_463", "TSG")
analyse_both("top_463", "oncogene")
analyse_both("top_463", "intogen")
analyse_both("top_463", "census")
analyse_both("top_463", "NCG")

print("---------------")

analyse_both("random_half_common", "CCGD")
analyse_both("random_half_common", "TSG")
analyse_both("random_half_common", "oncogene")
analyse_both("random_half_common", "intogen")
analyse_both("random_half_common", "census")
analyse_both("random_half_common", "NCG")

print("---------------")

overlap_matrix_random = []

print("Number of MFVS:", len(gene_lists))
for i in gene_lists:
    temp = []
    temp.append(len(list(j for j in i if j in database_dict["CCGD"])))
    temp.append(len(list(j for j in i if j in database_dict["census"])))
    temp.append(len(list(j for j in i if j in database_dict["intogen"])))
    temp.append(len(list(j for j in i if j in database_dict["oncogene"])))
    temp.append(len(list(j for j in i if j in database_dict["TSG"])))
    temp.append(len(list(j for j in i if j in database_dict["NCG"])))
    overlap_matrix_random.append(temp)

random_mean = np.mean(overlap_matrix_random, axis=0)
print("Average random MFVS in all genes:")
print("CCGD", random_mean[0])
print("census", random_mean[1])
print("intogen", random_mean[2])
print("oncogene", random_mean[3])
print("TSG", random_mean[4])
print("NCG", random_mean[5])
# print(overlap_matrix_random[0:4], overlap_matrix_random[-1])

def open_FVS_file(file_name, result_name):
    with open("FVS_result/{}.txt".format(file_name)) as f2:
        database_dict[result_name] = []
        for n, i in enumerate(f2):
            if n == 0:
                continue
            i = i.split()[0]
            database_dict[result_name].append(i)
    print("Size of {}:".format(result_name), len(database_dict[result_name]))


open_FVS_file("ILP_result_no_weight", "MFVS")
open_FVS_file("ILP_result_with_weight", "WMFVS")
open_FVS_file("ILP_result_with_weight_no_limitation", "WFVS")

print("---------------")

overlap_matrix_FVS = []
for i in [database_dict["MFVS"], database_dict["WMFVS"], database_dict["WFVS"]]:
    temp = []
    temp.append(len(list(j for j in i if j in database_dict["CCGD"])))
    temp.append(len(list(j for j in i if j in database_dict["census"])))
    temp.append(len(list(j for j in i if j in database_dict["intogen"])))
    temp.append(len(list(j for j in i if j in database_dict["oncogene"])))
    temp.append(len(list(j for j in i if j in database_dict["TSG"])))
    temp.append(len(list(j for j in i if j in database_dict["NCG"])))
    overlap_matrix_FVS.append(temp)

print(overlap_matrix_FVS)

overlap_matrix_random = np.array(overlap_matrix_random)
overlap_matrix_FVS = np.array(overlap_matrix_FVS)

CCGD_column = overlap_matrix_random[:, 0]
census_column = overlap_matrix_random[:, 1]
intogen_column = overlap_matrix_random[:, 2]
oncogene_column = overlap_matrix_random[:, 3]
TSG_column = overlap_matrix_random[:, 4]
NCG_column = overlap_matrix_random[:, 5]


def print_p_value_WMFVS(database_name, database_column, overlap_matrix_FVS_array):
    print("P-value of WMFVS in {}:".format(database_name))
    print("Max:", max(database_column), "Min:", min(database_column))
    print(overlap_matrix_FVS_array)
    temp_p_value = len(list(i for i in database_column if i >= overlap_matrix_FVS_array)) / len(database_column)
    print("P-value of {}: {}".format(database_name, temp_p_value))
    print("\n")

def print_rank_WMFVS(database_name, database_column, overlap_matrix_FVS_array):
    print("Rank of WMFVS in {}:".format(database_name))
    print("Max:", max(database_column), "Min:", min(database_column))
    print(overlap_matrix_FVS_array)
    temp_rank = len(list(i for i in database_column if i <= overlap_matrix_FVS_array)) / len(database_column)
    print("Rank of {}: {}".format(database_name, temp_rank))
    print("\n")

print_p_value_WMFVS("CCGD", CCGD_column, overlap_matrix_FVS[1][0])
print_p_value_WMFVS("census", census_column, overlap_matrix_FVS[1][1])
print_p_value_WMFVS("intogen", intogen_column, overlap_matrix_FVS[1][2])
print_p_value_WMFVS("oncogene", oncogene_column, overlap_matrix_FVS[1][3])
print_p_value_WMFVS("TSG", TSG_column, overlap_matrix_FVS[1][4])
print_p_value_WMFVS("NCG", NCG_column, overlap_matrix_FVS[1][5])

print("---------------")

print_rank_WMFVS("CCGD", CCGD_column, overlap_matrix_FVS[1][0])
print_rank_WMFVS("census", census_column, overlap_matrix_FVS[1][1])
print_rank_WMFVS("intogen", intogen_column, overlap_matrix_FVS[1][2])
print_rank_WMFVS("oncogene", oncogene_column, overlap_matrix_FVS[1][3])
print_rank_WMFVS("TSG", TSG_column, overlap_matrix_FVS[1][4])
print_rank_WMFVS("NCG", NCG_column, overlap_matrix_FVS[1][5])

print("---------------")

WFVS_resize = overlap_matrix_FVS[2, :]
WFVS_resize = WFVS_resize / 528 * 463
print("WFVS resized:", WFVS_resize)

def print_p_value_WFVS(database_name, database_column, overlap_matrix_FVS_array):
    print("P-value of WFVS in {}:".format(database_name))
    print("Max:", max(database_column), "Min:", min(database_column))
    print(overlap_matrix_FVS_array)
    temp_p_value = len(list(i for i in database_column if i >= overlap_matrix_FVS_array)) / len(database_column)
    print("P-value of {}: {}".format(database_name, temp_p_value))
    print("\n")

def print_rank_WFVS(database_name, database_column, overlap_matrix_FVS_array):
    print("Rank of WFVS in {}:".format(database_name))
    print("Max:", max(database_column), "Min:", min(database_column))
    print(overlap_matrix_FVS_array)
    temp_rank = len(list(i for i in database_column if i <= overlap_matrix_FVS_array)) / len(database_column)
    print("Rank of {}: {}".format(database_name, temp_rank))
    print("\n")

print_p_value_WFVS("CCGD", CCGD_column, WFVS_resize[0])
print_p_value_WFVS("census", census_column, WFVS_resize[1])
print_p_value_WFVS("intogen", intogen_column, WFVS_resize[2])
print_p_value_WFVS("oncogene", oncogene_column, WFVS_resize[3])
print_p_value_WFVS("TSG", TSG_column, WFVS_resize[4])
print_p_value_WFVS("NCG", NCG_column, WFVS_resize[5])

print("---------------")

print_rank_WFVS("CCGD", CCGD_column, WFVS_resize[0])
print_rank_WFVS("census", census_column, WFVS_resize[1])
print_rank_WFVS("intogen", intogen_column, WFVS_resize[2])
print_rank_WFVS("oncogene", oncogene_column, WFVS_resize[3])
print_rank_WFVS("TSG", TSG_column, WFVS_resize[4])
print_rank_WFVS("NCG", NCG_column, WFVS_resize[5])



print("---------------")

count = [0, 0, 0, 0, 0, 0]
for i in range(100):
    random.shuffle(all_gene_list)
    count[0] += len(list(i for i in database_dict["CCGD"] if i in all_gene_list[:463]))
    count[1] += len(list(i for i in database_dict["intogen"] if i in all_gene_list[:463]))
    count[2] += len(list(i for i in database_dict["census"] if i in all_gene_list[:463]))
    count[3] += len(list(i for i in database_dict["oncogene"] if i in all_gene_list[:463]))
    count[4] += len(list(i for i in database_dict["TSG"] if i in all_gene_list[:463]))
    count[5] += len(list(i for i in database_dict["NCG"] if i in all_gene_list[:463]))

count = list(i/100 for i in count)

print("Both in CCGD and random_463:")
print("Size:", count[0], "\n")
print("Precision:", count[0] / 463)
print("\n")

print("Both in intogen and random_463:")
print("Size:", count[1], "\n")
print("Precision:", count[1] / 463)
print("\n")

print("Both in census and random_463:")
print("Size:", count[2], "\n")
print("Precision:", count[2] / 463)
print("\n")

print("Both in oncogene and random_463:")
print("Size:", count[3], "\n")
print("Precision:", count[3] / 463)
print("\n")

print("Both in TSG and random_463:")
print("Size:", count[4], "\n")
print("Precision:", count[4] / 463)
print("\n")

print("Both in NCG and random_463:")
print("Size:", count[5], "\n")
print("Precision:", count[5] / 463)
print("\n")
