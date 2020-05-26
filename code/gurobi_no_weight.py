from gurobipy import *
import time

class Node:
    def __init__(self):
        self.in_list = []
        self.out_list = []
        self.in_degree = 0
        self.out_degree = 0
        self.label = ""
        self.name = ""

edge_list = []
node_list = []
node_list_count = 0
node_label_dict = {}

# with open("test_graph1.txt") as f1:
with open("result/compressed_graph.txt") as f1:
# with open("result/compressed_graph_no_ignore.txt") as f1:
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
                node_label_dict[node2_name] = node_list_count
                node_list_count += 1
                node_list.append(temp)
            node_list[node_label_dict[node1_name]].out_list.append(node2_name)
            node_list[node_label_dict[node1_name]].out_degree += 1
            node_list[node_label_dict[node2_name]].in_list.append(node1_name)
            node_list[node_label_dict[node2_name]].in_degree += 1
            edge_list.append((node_label_dict[node1_name], node_label_dict[node2_name]))

print("Read file finished, start ILP...")

time_start = time.clock()
model_1 = Model("test1")
variables_x = []
variables_w = []

for i in range(node_list_count):
    variables_x.append(model_1.addVar(vtype=GRB.BINARY, name="x{}".format(i)))
    variables_w.append(model_1.addVar(vtype=GRB.INTEGER, name='w{}'.format(i), lb=0, ub=node_list_count - 1))

model_1.setObjective(quicksum(variables_x), GRB.MINIMIZE)
for n, i in enumerate(edge_list):
    model_1.addConstr(variables_w[i[0]] - variables_w[i[1]] + node_list_count * variables_x[i[0]] >= 1, "c{}".format(i))
model_1.optimize()

time_end = time.clock()


#####################

result_list = []
with open("result/temp_MFVS.txt") as f2:
    for n, i in enumerate(f2):
        if n == 0:
            continue
        i = i.split()
        temp = Node()
        temp.name = i[0]
        temp.label = i[1]
        result_list.append(temp)

for v in model_1.getVars():
    if v.varName[0] == "x":
        if v.x == 1:
            index = int(v.varName[1:])
            result_list.append(node_list[index])

print('\n----Obj: {}----'.format(model_1.objVal))
print("Time used: {}".format(time_end-time_start))

with open("result/ILP_result_no_weight.txt", "w") as f2:
    f2.write("Number of MFVS nodes: {}\n".format(len(result_list)))
    for i in result_list:
        f2.write("{}\t{}\n".format(i.name, i.label))

print("end")
