import os
import numpy as np
import random

database_dict = {}

CCGD_genes = []
census_genes = []
intogen_genes = []
oncogene_genes = []
TSG_genes = []
NCG_genes = []
CGN_genes = []

gene_dict = {}
file_name_list = []
gene_lists = []
random_common_genes = []
random_half_common_genes = []
path = "random_result"

gurobi_MFVS_genes = []
gurobi_WMFVS_genes = []
gurobi_WFVS_genes = []

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

print("Got {} random files.".format(len(file_name_list)))
print("Got {} common genes.".format(len(random_common_genes)))
print("Got {} genes appearing in 49.5% MFVSs".format(len(random_half_common_genes)))


def open_gene_file(database_name):
    database_dict[database_name] = []
    with open("database/{}_genes.txt".format(database_name)) as f1:
        for n, i in enumerate(f1):
            i = i.split()
            database_dict[database_name].append(i[0])
    print("Size of {}_genes:".format(database_name), len(database_dict[database_name]))


open_gene_file("CCGD")

with open("database/CCGD_genes.txt") as f1:
    for n, i in enumerate(f1):
        i = i.split()
        CCGD_genes.append(i[0])
print("Size of CCGD_genes:", len(CCGD_genes))
''''''

with open("database/census_genes.txt") as f1:
    for n, i in enumerate(f1):
        i = i.split()
        census_genes.append(i[0])
print("Size of census_genes:", len(census_genes))

with open("database/intogen_genes.txt") as f1:
    for n, i in enumerate(f1):
        i = i.split()
        intogen_genes.append(i[0])
print("Size of intogen_genes:", len(intogen_genes))

with open("database/oncogene_genes.txt") as f1:
    for n, i in enumerate(f1):
        i = i.split()
        oncogene_genes.append(i[0])
print("Size of oncogene_genes:", len(oncogene_genes))

with open("database/TSG_genes.txt") as f1:
    for n, i in enumerate(f1):
        i = i.split()
        TSG_genes.append(i[0])
print("Size of TSG_genes:", len(TSG_genes))

with open("database/NCG_genes.txt") as f1:
    for n, i in enumerate(f1):
        i = i.split()
        NCG_genes.append(i[0])
print("Size of NCG_genes:", len(NCG_genes))

with open("gene_weight.txt") as f1:
    for n, i in enumerate(f1):
        i = i.split()
        gene_weight.append((i[0], float(i[-1])))

gene_weight = sorted(gene_weight, key=lambda x: x[-1])
top_463_genes = list(i[0] for i in gene_weight[-463:])

print("---------------")

print("Both in CCGD_genes and random_common_genes:")
temp = list(i for i in CCGD_genes if i in random_common_genes)
# print(sorted(temp))
print("Size:", len(temp), "\n")

print("Both in TSG_genes and random_common_genes:")
temp = list(i for i in TSG_genes if i in random_common_genes)
# print(sorted(temp))
print("Size:", len(temp), "\n")

print("Both in oncogene_genes and random_common_genes:")
temp = list(i for i in oncogene_genes if i in random_common_genes)
# print(sorted(temp))
print("Size:", len(temp), "\n")

print("Both in intogen_genes and random_common_genes:")
temp = list(i for i in intogen_genes if i in random_common_genes)
# print(sorted(temp))
print("Size:", len(temp), "\n")

print("Both in census and random_common_genes:")
temp = list(i for i in census_genes if i in random_common_genes)
# print(sorted(temp))
print("Size:", len(temp), "\n")

print("Both in NCG and random_common_genes:")
temp = list(i for i in NCG_genes if i in random_common_genes)
# print(sorted(temp))
print("Size:", len(temp), "\n")

print("---------------")

print("Both in CCGD_genes and random_half_common_genes:")
temp = list(i for i in CCGD_genes if i in random_half_common_genes)
# print(sorted(temp))
print("Size:", len(temp), "\n")

print("Both in TSG_genes and random_half_common_genes:")
temp = list(i for i in TSG_genes if i in random_half_common_genes)
# print(sorted(temp))
print("Size:", len(temp), "\n")

print("Both in oncogene_genes and random_half_common_genes:")
temp = list(i for i in oncogene_genes if i in random_half_common_genes)
# print(sorted(temp))
print("Size:", len(temp), "\n")

print("Both in intogen_genes and random_half_common_genes:")
temp = list(i for i in intogen_genes if i in random_half_common_genes)
# print(sorted(temp))
print("Size:", len(temp), "\n")

print("Both in census and random_half_common_genes:")
temp = list(i for i in census_genes if i in random_half_common_genes)
# print(sorted(temp))
print("Size:", len(temp), "\n")

print("Both in NCG and random_half_common_genes:")
temp = list(i for i in NCG_genes if i in random_half_common_genes)
# print(sorted(temp))
print("Size:", len(temp), "\n")

print("---------------")

overlap_matrix_random = []

print("Number of MFVS:", len(gene_lists))
for i in gene_lists:
    temp = []
    temp.append(len(list(j for j in i if j in CCGD_genes)))
    temp.append(len(list(j for j in i if j in census_genes)))
    temp.append(len(list(j for j in i if j in intogen_genes)))
    temp.append(len(list(j for j in i if j in oncogene_genes)))
    temp.append(len(list(j for j in i if j in TSG_genes)))
    temp.append(len(list(j for j in i if j in NCG_genes)))
    overlap_matrix_random.append(temp)

# print(overlap_matrix_random[0:4], overlap_matrix_random[-1])

with open("FVS_result/ILP_result_no_weight.txt") as f2:
    for n, i in enumerate(f2):
        if n == 0:
            continue
        i = i.split()[0]
        gurobi_MFVS_genes.append(i)
print("Size of gurobi_no_weight_genes:", len(gurobi_MFVS_genes))

with open("FVS_result/ILP_result_with_weight.txt") as f2:
    for n, i in enumerate(f2):
        if n == 0:
            continue
        i = i.split()[0]
        gurobi_WMFVS_genes.append(i)
print("Size of gurobi_with_weight_genes:", len(gurobi_WMFVS_genes))

with open("FVS_result/ILP_result_with_weight_no_limitation.txt") as f2:
    for n, i in enumerate(f2):
        if n == 0:
            continue
        i = i.split()[0]
        gurobi_WFVS_genes.append(i)
print("Size of gurobi_with_weight_no_limitation_genes:", len(gurobi_WFVS_genes))

print("---------------")

overlap_matrix_FVS = []
for i in [gurobi_MFVS_genes, gurobi_WMFVS_genes, gurobi_WFVS_genes]:
    temp = []
    temp.append(len(list(j for j in i if j in CCGD_genes)))
    temp.append(len(list(j for j in i if j in census_genes)))
    temp.append(len(list(j for j in i if j in intogen_genes)))
    temp.append(len(list(j for j in i if j in oncogene_genes)))
    temp.append(len(list(j for j in i if j in TSG_genes)))
    temp.append(len(list(j for j in i if j in NCG_genes)))
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


print("Rank of WMFVS in CCGD:")
# print(CCGD_column)
print("Max:", max(CCGD_column), "Min:", min(CCGD_column))
print(overlap_matrix_FVS[1][0])
# print(list(i for i in CCGD_column if i <= overlap_matrix_FVS[1][0]))
CCGD_WMFVS_rate = len(list(i for i in CCGD_column if i <= overlap_matrix_FVS[1][0])) / len(CCGD_column)
print(CCGD_WMFVS_rate)
print("\n")

print("Rank of WMFVS in census:")
# print(census_column)
print(overlap_matrix_FVS[1][1])
print("Max:", max(census_column), "Min:", min(census_column))
# print(list(i for i in census_column if i <= overlap_matrix_FVS[1][1]))
census_WMFVS_rate = len(list(i for i in census_column if i <= overlap_matrix_FVS[1][1])) / len(census_column)
print(census_WMFVS_rate)
print("\n")

print("Rank of WMFVS in intogen:")
# print(intogen_column)
print(overlap_matrix_FVS[1][2])
print("Max:", max(intogen_column), "Min:", min(intogen_column))
# print(list(i for i in intogen_column if i <= overlap_matrix_FVS[1][2]))
intogen_WMFVS_rate = len(list(i for i in intogen_column if i <= overlap_matrix_FVS[1][2])) / len(intogen_column)
print(intogen_WMFVS_rate)
print("\n")

print("Rank of WMFVS in oncogene:")
# print(oncogene_column)
print(overlap_matrix_FVS[1][3])
print("Max:", max(oncogene_column), "Min:", min(oncogene_column))
print(sorted(list(i for i in oncogene_column if i <= overlap_matrix_FVS[1][3])))
print(len(sorted(list(i for i in oncogene_column if i <= overlap_matrix_FVS[1][3]))))
oncogene_WMFVS_rate = len(list(i for i in oncogene_column if i <= overlap_matrix_FVS[1][3])) / len(oncogene_column)
print(oncogene_WMFVS_rate)
print("\n")

print("Rank of WMFVS in TSG:")
# print(TSG_column)
print(overlap_matrix_FVS[1][4])
print("Max:", max(TSG_column), "Min:", min(TSG_column))
print(sorted(list(i for i in TSG_column if i <= overlap_matrix_FVS[1][4])))
print(len(sorted(list(i for i in TSG_column if i <= overlap_matrix_FVS[1][4]))))
TSG_WMFVS_rate = len(list(i for i in TSG_column if i <= overlap_matrix_FVS[1][4])) / len(TSG_column)
print(TSG_WMFVS_rate)
print("\n")

print("Rank of WMFVS in NCG:")
# print(NCG_column)
print(overlap_matrix_FVS[1][5])
print("Max:", max(NCG_column), "Min:", min(NCG_column))
print(sorted(list(i for i in NCG_column if i <= overlap_matrix_FVS[1][5])))
print(len(sorted(list(i for i in NCG_column if i <= overlap_matrix_FVS[1][5]))))
NCG_WMFVS_rate = len(list(i for i in NCG_column if i <= overlap_matrix_FVS[1][5])) / len(NCG_column)
print(NCG_WMFVS_rate)
print("\n")

print("---------------")
WFVS_resize = overlap_matrix_FVS[2, :]
WFVS_resize = WFVS_resize / 528 * 463
print("WFVS resized:", WFVS_resize)

print("Rank of WFVS in CCGD:")
# print(CCGD_column)
print("Max:", max(CCGD_column), "Min:", min(CCGD_column))
print(WFVS_resize[0])
# print(list(i for i in CCGD_column if i <= overlap_matrix_FVS[1][0]))
CCGD_WFVS_rate = len(list(i for i in CCGD_column if i <= WFVS_resize[0])) / len(CCGD_column)
print(CCGD_WFVS_rate)
print("\n")

print("Rank of WFVS in census:")
# print(census_column)
print(WFVS_resize[1])
print("Max:", max(census_column), "Min:", min(census_column))
# print(list(i for i in census_column if i <= overlap_matrix_FVS[1][1]))
census_WFVS_rate = len(list(i for i in census_column if i <= WFVS_resize[1])) / len(census_column)
print(census_WFVS_rate)
print("\n")

print("Rank of WFVS in intogen:")
# print(intogen_column)
print(WFVS_resize[2])
print("Max:", max(intogen_column), "Min:", min(intogen_column))
# print(list(i for i in intogen_column if i <= overlap_matrix_FVS[1][2]))
intogen_WFVS_rate = len(list(i for i in intogen_column if i <= WFVS_resize[2])) / len(intogen_column)
print(intogen_WFVS_rate)
print("\n")

print("Rank of WFVS in oncogene:")
# print(oncogene_column)
print(WFVS_resize[3])
print("Max:", max(oncogene_column), "Min:", min(oncogene_column))
print(sorted(list(i for i in oncogene_column if i <= WFVS_resize[3])))
print(len(sorted(list(i for i in oncogene_column if i <= WFVS_resize[3]))))
oncogene_WFVS_rate = len(list(i for i in oncogene_column if i <= WFVS_resize[3])) / len(oncogene_column)
print(oncogene_WFVS_rate)
print("\n")

print("Rank of WFVS in TSG:")
# print(TSG_column)
print(WFVS_resize[4])
print("Max:", max(TSG_column), "Min:", min(TSG_column))
print(sorted(list(i for i in TSG_column if i <= WFVS_resize[4])))
print(len(sorted(list(i for i in TSG_column if i <= WFVS_resize[4]))))
TSG_WFVS_rate = len(list(i for i in TSG_column if i <= WFVS_resize[4])) / len(TSG_column)
print(TSG_WFVS_rate)
print("\n")

print("Rank of WFVS in NCG:")
# print(NCG_column)
print(WFVS_resize[5])
print("Max:", max(NCG_column), "Min:", min(NCG_column))
print(sorted(list(i for i in NCG_column if i <= WFVS_resize[5])))
print(len(sorted(list(i for i in NCG_column if i <= WFVS_resize[5]))))
NCG_WFVS_rate = len(list(i for i in NCG_column if i <= WFVS_resize[5])) / len(NCG_column)
print(NCG_WFVS_rate)
print("\n")

print("---------------")

print("Both in CCGD and top_463:")
temp = list(i for i in CCGD_genes if i in top_463_genes)
print(sorted(temp))
print("Size:", len(temp), "\n")
print("\n")

print("Both in intogen and top_463:")
temp = list(i for i in intogen_genes if i in top_463_genes)
print(sorted(temp))
print("Size:", len(temp), "\n")
print("\n")

print("Both in census and top_463:")
temp = list(i for i in census_genes if i in top_463_genes)
print(sorted(temp))
print("Size:", len(temp), "\n")
print("\n")

print("Both in oncogene and top_463:")
temp = list(i for i in oncogene_genes if i in top_463_genes)
print(sorted(temp))
print("Size:", len(temp), "\n")
print("\n")

print("Both in TSG and top_463:")
temp = list(i for i in TSG_genes if i in top_463_genes)
print(sorted(temp))
print("Size:", len(temp), "\n")
print("\n")

print("Both in NCG and top_463:")
temp = list(i for i in NCG_genes if i in top_463_genes)
print(sorted(temp))
print("Size:", len(temp), "\n")
print("\n")

print("---------------")
all_gene_list = list(i[0] for i in gene_weight)

count = [0, 0, 0, 0, 0, 0]
for i in range(100):
    random.shuffle(all_gene_list)
    count[0] += len(list(i for i in CCGD_genes if i in all_gene_list[:463]))
    count[1] += len(list(i for i in intogen_genes if i in all_gene_list[:463]))
    count[2] += len(list(i for i in census_genes if i in all_gene_list[:463]))
    count[3] += len(list(i for i in oncogene_genes if i in all_gene_list[:463]))
    count[4] += len(list(i for i in TSG_genes if i in all_gene_list[:463]))
    count[5] += len(list(i for i in NCG_genes if i in all_gene_list[:463]))

count = list(i/100 for i in count)

print("Both in CCGD and random_463:")
print("Size:", count[0], "\n")
print("\n")

print("Both in intogen and top_463:")
print("Size:", count[1], "\n")
print("\n")

print("Both in census and top_463:")
print("Size:", count[2], "\n")
print("\n")

print("Both in oncogene and top_463:")
print("Size:", count[3], "\n")
print("\n")

print("Both in TSG and top_463:")
print("Size:", count[4], "\n")
print("\n")

print("Both in NCG and top_463:")
print("Size:", count[5], "\n")
print("\n")
