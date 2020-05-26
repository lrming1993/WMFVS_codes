genes = []
with open("result/gene_weight.txt") as f1:
    for i in f1:
        i = i.split()
        genes.append((i[0], float(i[-1])))

genes = sorted(genes, key=lambda x: x[1], reverse=True)
print(genes)

with open("result/gene_weight_sorted.txt", "w") as f2:
    f2.write("GENE\tabs_logFC\n")
    for i in genes:
        f2.write("{}\t{}\n".format(i[0], i[1]))

TSG_genes = []

with open("resource/cancer_genes/Human_TSGs.txt") as f1:
    for n, i in enumerate(f1):
        if n == 0:
            continue
        i = i.split()
        TSG_genes.append(i[1])

print(TSG_genes)

with open("result/phenotype_TSG.cls", "w") as f2:
    f2.write("{} {} {}\n".format(len(genes), 2, 1))
    f2.write("# True False\n")
    for n, i in enumerate(genes):
        f2.write("{}{}".format("1" if i[0] in TSG_genes else "0", " " if n < len(genes) - 1 else ""))
