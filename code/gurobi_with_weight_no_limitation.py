from gurobipy import *
import time

penalty_bound = 65536

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
        self.penalty = penalty_bound

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
temp_list = []
with open("result/temp_MFVS_no_ignore.txt") as f2:
    for n, i in enumerate(f2):
        if n == 0:
            continue
        i = i.split()
        temp = Node()
        temp.name = i[0]
        temp.label = i[1]
        temp_list.append(temp)

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

with open("result/compressed_graph_no_ignore.txt") as f1:
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
        i.penalty = 1 / i.weight if abs(i.weight) >= 1 / penalty_bound else penalty_bound

print("Read file finished, start ILP...")


time_start = time.clock()
model_1 = Model("test1")
variables_x = []
variables_w = []

for i in range(node_list_count):
    variables_x.append(model_1.addVar(vtype=GRB.BINARY, name="x{}".format(i)))
    variables_w.append(model_1.addVar(vtype=GRB.INTEGER, name='w{}'.format(i), lb=0, ub=node_list_count - 1))

model_1.setObjective(quicksum(i.penalty * variables_x[i.index] for i in node_list), GRB.MINIMIZE)

for n, i in enumerate(edge_list):
    model_1.addConstr(variables_w[i[0]] - variables_w[i[1]] + node_list_count * variables_x[i[0]] >= 1, "c{}".format(i))

model_1.optimize()

time_end = time.clock()


#####################

result_list = []
for v in model_1.getVars():
    if v.varName[0] == "x":
        if v.x == 1:
            index = int(v.varName[1:])
            result_list.append(node_list[index])

result_list += temp_list

print('\n----Obj: {}----'.format(model_1.objVal))
print("Time used: {}".format(time_end-time_start))

with open("result/ILP_result_with_weight_no_limitation.txt", "w") as f2:
    f2.write("Number of FVS nodes: {}\n".format(len(result_list)))
    for i in result_list:
        f2.write("{}\t{}\t{}\n".format(i.name, i.label, i.weight))

print("end")
