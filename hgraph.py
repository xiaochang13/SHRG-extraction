#!/usr/bin/env python3
import gflags
FLAGS = gflags.FLAGS

class Hyperedge(object):
    """
    A class for hyperedge. Each hyperedge consists of a tuple
    of ordered vertices. Each hyperedge also has an edge label.
    """
    def __init__(self, vertices_, label_, weight_=1.0):
        self.vertices = vertices_
        self.label = label_
        self.weight = weight_

    def degree(self):
        return len(self.vertices)

    def __str__(self):
        return "(%s:%s)" % (",".join(self.vertices), self.label)

class Hypergraph(object):
    def __init__(self):
        self.vertices = []  # A list of variable names for all vertices.
        self.vertex_idx = {}
        self.edges = []  # A list of hyperedges, ordered according to vertex order.
        self.ext_list = []  # An ordered list of vertices.
        self.lhs = None

    def init_vertices(self, bag):
        bag_list = sorted([id for id in bag])
        for vertex in bag_list:
            assert vertex not in self.vertex_idx
            curr_lab = "n%d" % len(self.vertices)
            self.vertex_idx[vertex] = curr_lab
            self.vertices.append(curr_lab)

    def setExtList(self, ext_list_):
        self.ext_list = ext_list_

    def addEdge(self, edge):
        self.edges.append(edge)

    def __str__(self):
        edge_repr = ";".join([str(edge) for edge in self.edges])
        ext_repr = ",".join(self.ext_list)
        return "X%d -> %s | %s" % (len(self.ext_list), edge_repr, ext_repr)

