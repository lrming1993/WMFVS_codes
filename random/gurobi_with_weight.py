from gurobipy import *
import time
import random
import os

MFVS_num = 463

class Node:
    def __init__(self):
        self.in_list = []
        self.out_list = []
        self.in_degree = 0
        self.out_degree = 0
        self.label = "-"
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
gene_symbol = []
gene_FC = []
gene_symbol_dict = {}

with open("gene_weight.txt") as f1:
    for n, i in enumerate(f1):
        i = i.split()
        gene_symbol.append(i[0])
        gene_FC.append(float(i[2]))
        gene_symbol_dict[i[0]] = n

print(len(gene_symbol))
print(len(gene_FC))


with open("compressed_graph_no_ignore.txt") as f1:
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

print(len(node_list))

for _ in range(1000):
    random.shuffle(gene_FC)
    gene_list = []

    for n, i in enumerate(gene_symbol):
        temp = Gene()
        temp.symbol = i
        temp.FC = gene_FC[n]
        temp.index = n
        gene_list.append(temp)

    for i in node_list:
        if i.name in gene_symbol_dict:
            i.weight = abs(gene_list[gene_symbol_dict[i.name]].FC)

    print("Read file finished, start ILP...")

    time_start = time.clock()
    model_1 = Model("model_{}".format(_))
    variables_x = []
    variables_w = []

    for i in range(node_list_count):
        variables_x.append(model_1.addVar(vtype=GRB.BINARY, name="x{}".format(i)))
        variables_w.append(model_1.addVar(vtype=GRB.INTEGER, name='w{}'.format(i), lb=0, ub=node_list_count - 1))

    model_1.setObjective(quicksum(i.weight * variables_x[i.index] for i in node_list), GRB.MAXIMIZE)

    # test = "c{}".format(edge_list[0])
    for n, i in enumerate(edge_list):
        model_1.addConstr(variables_w[i[0]] - variables_w[i[1]] + node_list_count * variables_x[i[0]] >= 1, "c{}".format(n))
    model_1.addConstr(quicksum(variables_x) >= MFVS_num)
    model_1.addConstr(quicksum(variables_x) <= MFVS_num)

    model_1.optimize()

    time_end = time.clock()


    #####################

    result_list = []
    for v in model_1.getVars():
        if v.varName[0] == "x":
            if v.x == 1:
                index = int(v.varName[1:])
                result_list.append(node_list[index])
    print('\n----Obj: {}----'.format(model_1.objVal))
    print("Time used: {}".format(time_end-time_start))
    print("Size of result:", len(result_list))
    if len(result_list) != MFVS_num:
        temp_path = "result_0320/fault_set/{}".format(int(time.time()))
        os.mkdir(temp_path)
        with open("{}/gene_weight_fault_{}.txt".format(temp_path, int(time.time())), "w") as f3:
            for n, i in enumerate(gene_symbol):
                f3.write("{}\t{}\n".format(i, gene_FC[n]))
        with open("{}/ILP_result_with_random_weight_{}_{}.txt".format(temp_path, int(time.time()), _), "w") as f2:
            f2.write("Number of MFVS nodes: {}\n".format(len(result_list)))
            for i in result_list:
                f2.write("{}\t{}\t{}\n".format(i.name, i.label, i.weight))
    else:
        with open("result_0320/result/ILP_result_with_random_weight_{}_{}.txt".format(int(time.time()), _), "w") as f2:
            f2.write("Number of MFVS nodes: {}\n".format(len(result_list)))
            for i in result_list:
                f2.write("{}\t{}\t{}\n".format(i.name, i.label, i.weight))

    print("end")
