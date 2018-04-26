import sys, os
from dependency import DependencyTree
from collections import defaultdict
NULL = "-NULL-"
UNKNOWN = "-UNK-"

def readToks(path):
    tok_seqs = []
    with open(path, "r") as tok_f:
        for line in tok_f:
            line = line.strip()
            if not line:
                break
            toks = line.split()
            tok_seqs.append(toks)
    return tok_seqs

def loadDependency(path):
    dep_trees = []
    all_toks = []
    n_fields = -1
    num_vertices = 0
    edge_map = defaultdict(set)
    edge_label_map = {}
    toks = []
    with open(path, "r") as dep_f:
        for line in dep_f:
            splits = line.strip().split("\t")
            if len(splits) > 2:
                if n_fields == -1:
                    n_fields = len(splits)
                assert n_fields == len(splits)
            if len(splits) != n_fields:
                dep_trees.append((num_vertices, edge_map, edge_label_map))
                all_toks.append(toks)
                edge_map = defaultdict(set)
                edge_label_map = {}
                toks = []
            else:
                dep_label = splits[7]
                tok_idx = int(splits[0]) - 1
                toks.append(splits[1])
                head_idx = int(splits[6]) - 1 # root becomes -1.
                num_vertices = tok_idx + 1

                if head_idx != -1:
                    edge_map[head_idx].add(tok_idx)
                    edge_map[tok_idx].add(head_idx)
                    if tok_idx < head_idx: #
                        arc_lab = "L-" + dep_label
                        edge_label_map[(tok_idx, head_idx)] = arc_lab
                    else:
                        arc_lab = "R-" + dep_label
                        edge_label_map[(head_idx, tok_idx)] = arc_lab

        dep_f.close()
    # for triple in dep_trees:
    #     print triple[0], triple[1], triple[2]
    return all_toks, dep_trees

# loadDependency(sys.argv[1])
