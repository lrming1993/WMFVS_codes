import time
import copy

class Node:
    def __init__(self):
        self.in_list = []
        self.out_list = []
        self.in_degree = 0
        self.out_degree = 0
        self.label = ""
        self.name = ""
        self.flag = 1
        self.index = -1
        self.group = -1

edge_list = []
node_list = []
# compressed_node_list = []
# compressed_edge_list = []
node_list_count = 0
node_name_dict = {}
# compressed_node_name_dict = {}
# compress_mapping_dict = {}  # 记录压缩后节点index的变化，key为压缩前index，value为压缩后index

result_list = []

gene_weight_dict = {}

with open("result/gene_weight.txt") as f1:
    for i in f1:
        i = i.split()
        gene_weight_dict[i[0]] = float(i[-1])

def edge_add(node1, node2):
    global edge_list
    global node_list
    node_list[node1.index].out_list.append(node2.name)
    node_list[node1.index].out_degree += 1
    node_list[node2.index].in_list.append(node1.name)
    node_list[node2.index].in_degree += 1
    edge_list.append((node1.index, node2.index))


def node_deletion(node_index):
    global node_list
    global edge_list
    node_list[node_index].flag = 0
    for i in node_list[node_index].in_list:
        node_list[node_name_dict[i]].out_list.remove(node_list[node_index].name)
        node_list[node_name_dict[i]].out_degree -= 1
    for i in node_list[node_index].out_list:
        node_list[node_name_dict[i]].in_list.remove(node_list[node_index].name)
        node_list[node_name_dict[i]].in_degree -= 1
    edge_list = list(i for i in edge_list if node_index not in i)


def node_ignore(node_index):
    global node_list
    global edge_list
    in_list = node_list[node_index].in_list
    out_list = node_list[node_index].out_list
    for i in in_list:
        for j in out_list:
            edge_add(node_list[node_name_dict[i]], node_list[node_name_dict[j]])
    node_deletion(node_index)


def edge_deletion(node1_index, node2_index):
    global edge_list
    global node_list
    edge_list.remove((node1_index, node2_index))
    node_list[node1_index].out_list.remove(node_list[node2_index].name)
    node_list[node1_index].out_degree -= 1
    node_list[node2_index].in_list.remove(node_list[node1_index].name)
    node_list[node2_index].in_degree -= 1


# parameters for Gabow algorithm
order = [-1] * len(node_list)
part = [-1] * len(node_list)
path = []
root = []
order_number = 0
part_number = 0


def gabow(node_index):
    global order
    global part
    global path
    global root
    global order_number
    global part_number
    order_number += 1
    order[node_index] = order_number
    path.append(node_index)
    root.append(node_index)
    for i in node_list[node_index].out_list:
        if order[node_name_dict[i]] == -1:
            gabow(node_name_dict[i])
        elif part[node_name_dict[i]] == -1:
            while order[root[-1]] > order[node_name_dict[i]]:
                root = root[:-1]
    if node_index == root[-1]:
        part_number += 1
        root = root[:-1]
        while True:
            top = path[-1]
            part[top] = part_number
            path = path[:-1]
            if top == node_index:
                break

def scc():
    operation_flag = 0
    global order
    global part
    global path
    global root
    global order_number
    global part_number
    global node_list
    global edge_list
    order = [-1] * len(node_list)
    part = [-1] * len(node_list)
    path = []
    root = []
    order_number = 0
    part_number = 0
    for n, i in enumerate(order):
        if i == -1:
            gabow(n)
    for i in range(len(node_list)):
        node_list[i].group = part[i]
    for i in node_list:
        for j in i.out_list:
            j = node_name_dict[j]
            if i.group != node_list[j].group:
                edge_deletion(i.index, j)
                operation_flag = 1
    print("SCC part number:", part_number)
    return operation_flag


def compress_graph():
    global result_list
    global node_name_dict
    global node_list
    global edge_list
    operation_flag = 1
    _operation_flag = 0
    while operation_flag == 1:
        operation_flag = 0
        for i in node_list:
            if i.flag == 0:
                continue
            if i.name in i.out_list:
                temp = Node()
                temp.name = i.name
                temp.label = i.label
                result_list.append(temp)
                node_deletion(i.index)
                operation_flag = 1
                _operation_flag = 1
                continue
            if i.in_degree == 0 or i.out_degree == 0:
                node_deletion(i.index)
                operation_flag = 1
                _operation_flag = 1
                continue
            if i.in_degree == 1:
                pre = i.in_list[0]
                if gene_weight_dict[pre] > gene_weight_dict[i.name]:
                    node_ignore(i.index)
                    operation_flag = 1
                    _operation_flag = 1
                    continue
            if i.out_degree == 1:
                suc = i.out_list[0]
                if gene_weight_dict[suc] > gene_weight_dict[i.name]:
                    node_ignore(i.index)
                    operation_flag = 1
                    _operation_flag = 1
                    continue
    if _operation_flag:
        compress_mapping_dict = {}
        compressed_node_name_dict = {}
        compressed_edge_list = []
        compressed_node_list = list(copy.deepcopy(i) for i in node_list if i.flag == 1)
        for n, i in enumerate(compressed_node_list):
            compress_mapping_dict[i.index] = n
            i.index = n
            compressed_node_name_dict[i.name] = n
        for i in compressed_node_list:
            for j in i.out_list:
                compressed_edge_list.append((i.index, compressed_node_list[compressed_node_name_dict[j]].index))
        node_list = compressed_node_list
        edge_list = compressed_edge_list
        node_name_dict = compressed_node_name_dict
    return _operation_flag


# with open("test_graph3.txt") as f1:
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
            if node1_name not in node_name_dict:
                temp = Node()
                temp.label = node1_label
                temp.name = node1_name
                temp.index = node_list_count
                node_name_dict[node1_name] = node_list_count
                node_list_count += 1
                node_list.append(temp)
            if node2_name not in node_name_dict:
                temp = Node()
                temp.label = node2_label
                temp.name = node2_name
                temp.index = node_list_count
                node_name_dict[node2_name] = node_list_count
                node_list_count += 1
                node_list.append(temp)
            node_list[node_name_dict[node1_name]].out_list.append(node2_name)
            node_list[node_name_dict[node1_name]].out_degree += 1
            node_list[node_name_dict[node2_name]].in_list.append(node1_name)
            node_list[node_name_dict[node2_name]].in_degree += 1
            edge_list.append((node_name_dict[node1_name], node_name_dict[node2_name]))

print("Read file finished, start compress...")

# node_list.sort(key=lambda x: int(x.label))
count = 1
while True:
    print("Compress procedure {}".format(count))
    flag1 = scc()
    flag2 = compress_graph()
    count += 1
    if not(flag1 or flag2):
        print("No change.\n")
        break
    else:
        print("Remained node size:", len(node_list))
        print("Remained edge size:", len(edge_list), "\n")

print("----Compress finished!----")
print("Remained node size:", len(node_list))
print("Remained edge size:", len(edge_list))
print("Get {} MFVS nodes while compressing:\n{}".format(len(result_list), [i.name for i in result_list]))

with open("result/compressed_graph_weight_ignore.txt", "w") as f3:
    f3.write("node1_name\tnode1_label\tnode2_name\tnode2_label\n")
    for i in edge_list:
        f3.write("{}\t{}\t{}\t{}\n".format(node_list[i[0]].name, node_list[i[0]].label, node_list[i[1]].name, node_list[i[1]].label))

with open("result/temp_MFVS_weight_ignore.txt", "w") as f4:
    f4.write("Number of nodes: {}\n".format(len(result_list)))
    for i in result_list:
        f4.write("{}\t{}\n".format(i.name, i.label))
