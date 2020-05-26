class Node:
    def __init__(self):
        self.in_list = []
        self.out_list = []
        self.in_degree = 0
        self.out_degree = 0
        self.label = ""
        self.name = ""
        self.weight = 0.0
        self.index = 0

edge_list = []
node_list = []
node_list_count = 0
node_label_dict = {}

class Gene():
    def __init__(self):
        self.symbol = ""
        self.id = ""
        self.FC = ""
        self.p_value = ""
        self.index = -1

gene_list = []
gene_list_count = 0
gene_symbol_dict = {}
lack_list = []

with open("BRCA_merge_output_v2_sample_sort_NtoP_log2transformed_limma.txt") as f1:
    for n, i in enumerate(f1):
        if n == 0:
            continue
        i = i.split()
        temp = Gene()
        temp.symbol = i[0].split("|")[0]
        temp.id = i[0].split("|")[1]
        temp.FC = float(i[1])
        temp.p_value = float(i[5])
        temp.index = gene_list_count
        gene_symbol_dict[temp.symbol] = temp.index
        gene_list_count += 1
        gene_list.append(temp)

with open("2001699_Tables_S1_S2_S6.txt") as f1:
    for n, li in enumerate(f1):
        if n == 0:
            continue
        else:
            li = li.split()
            node1_name = li[0]
            node1_label = li[1]
            node2_name = li[2]
            node2_label = li[3]
            if node1_name not in node_label_dict:
                temp = Node()
                temp.in_list = []
                temp.out_list = []
                temp.in_degree = 0
                temp.out_degree = 0
                temp.label = node1_label
                temp.name = node1_name
                temp.index = node_list_count
                node_label_dict[node1_name] = node_list_count
                node_list_count += 1
                node_list.append(temp)
            if node2_name not in node_label_dict:
                temp = Node()
                temp.in_list = []
                temp.out_list = []
                temp.in_degree = 0
                temp.out_degree = 0
                temp.label = node2_label
                temp.name = node2_name
                temp.index = node_list_count
                node_label_dict[node2_name] = node_list_count
                node_list_count += 1
                node_list.append(temp)
            node_list[node_label_dict[node1_name]].out_list.append(node2_name)
            node_list[node_label_dict[node1_name]].out_degree += 1
            node_list[node_label_dict[node2_name]].in_list.append(node1_name)
            node_list[node_label_dict[node2_name]].in_degree += 1
            edge_list.append((node_label_dict[node1_name], node_label_dict[node2_name]))

for i in node_list:
    if i.name in gene_symbol_dict:
        i.weight = abs(gene_list[gene_symbol_dict[i.name]].FC)

with open("result/gene_weight.txt", "w") as f1:
    for i in node_list:
        f1.write("{}\t{}\t{}\n".format(i.name, i.label, i.weight))
