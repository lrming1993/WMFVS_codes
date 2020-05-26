import numpy as np
import matplotlib.pyplot as plt

graph_genes = []
gene_weight = []

CCGD_genes = []
census_genes = []
intogen_genes = []
oncogene_genes = []
TSG_genes = []

with open("gene_weight.txt") as f1:
    for n, i in enumerate(f1):
        i = i.split()
        graph_genes.append(i[0])
        gene_weight.append((i[0], float(i[2])))
gene_weight.sort(key=lambda x: x[1])
print("Size of graph_genes:", len(graph_genes))

with open("database/CCGD_genes.txt") as f1:
    for n, i in enumerate(f1):
        i = i.split()
        CCGD_genes.append(i[0])
print("Size of CCGD_genes:", len(CCGD_genes))

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

weight_list = np.array(gene_weight)[:,1].astype(np.float32)
gene_list = np.array(gene_weight)[:,0]
weight_list = (weight_list-weight_list.min())/(weight_list.max()-weight_list.min())
gene_weight = list(zip(gene_list, weight_list))
print(gene_weight)

def calculate_cov(database, msg):
    matrix = []
    for n, i in enumerate(gene_weight):
        if i[0] in database:
            matrix.append((i[1], 1))
        else:
            matrix.append((i[1], 0))
    matrix = np.array(matrix).T
    print("Cov of weight and {}".format(msg))
    print(np.cov(matrix))
    print("---------------")

print("---------------")
calculate_cov(CCGD_genes, "CCGD")
calculate_cov(census_genes, "census")
calculate_cov(intogen_genes, "intogen")
calculate_cov(oncogene_genes, "oncogene")
calculate_cov(TSG_genes, "TSG")


color_list = []
for i in gene_weight:
    color_list.append((1, 0, 0, 1) if i[0] in TSG_genes else (0, 0, 1, 0.1))

plt.scatter(np.arange(len(weight_list)), [1]*len(weight_list), c=color_list, marker=".")

plt.show()
