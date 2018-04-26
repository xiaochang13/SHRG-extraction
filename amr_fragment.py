#!/usr/bin/env python2.7
import sys
import bitarray
import copy
from bitarray import bitarray
from util.cfg import NonterminalLabel
import gflags
from HRGSample import Sample
#Each fragment is simply two vectors of bits, for edges and boundary nodes
class AMRFragment(object):
    def __init__(self, n_edges, n_nodes, graph):
        self.edges = bitarray(n_edges)
        if self.edges.count() != 0:
            self.edges ^= self.edges
        assert self.edges.count() == 0, 'initialization nonzero'
        self.nodes = bitarray(n_nodes)
        if self.nodes.count() != 0:
            self.nodes ^= self.nodes
        assert self.nodes.count() == 0, 'initialization nonzero'
        self.graph = graph
        #self.root = -1
        self.roots = []
        #self.roots = set() #This is used to allow multiple root structure
        self.ext_set = set()
        self.ext_list = []
        self.outside_set = set()
        self.inside_set = set()
        self.ext_mapping = {}
        self.start = -1
        self.end = -1
        self.unaligned = []
        self.ext_label = {} #This attribute is used for the refinement of the external nodes
        self.category = None
        self.toks = []

    def edgeCount(self):
        return self.edges.count()

    def allEdges(self):
        return set([i for i in xrange(len(self.edges)) if self.edges[i] == 1])

    def setSpan(self, start, end):
        self.start = start
        self.end = end
        if not self.category:
            self.toks = self.graph.sent[start:end]

    def setRoot(self, node_num):
        self.set_node(node_num)
        #self.root = node_num
        self.roots.append(node_num)

    def setCategory(self, cate_):
        self.category = cate_
        self.toks = [cate_]

    def rootEdge(self):
        assert len(self.roots) == 1
        root_index = self.roots[0]
        root_node = self.graph.nodes[root_index]
        return root_node.c_edge

    def str_side(self):
        if self.start == -1:
            return ''
        return ' '.join(self.graph.sent[self.start:self.end])

    def str_list(self):
        if self.start == -1:
            return []

        return self.graph.sent[self.start:self.end]

    def extSuffix(self, first_index):
        #assert first_index in self.ext_set, self.ext_set
        amr = self.graph
        ret = ""
        other_indices = []
        if first_index is not None and first_index in self.ext_set:
            const_edge = amr.nodes[first_index].c_edge
            ret = '1' if self.edges[const_edge] == 1 else '0'

        for index in self.ext_list:
            if index == first_index:
                continue
            other_indices.append(index)
            const_edge = amr.nodes[index].c_edge
            if self.edges[const_edge] == 1:
                ret += '1'
            else:
                ret += '0'
        return ret, other_indices

    #Return all (index, rel, index) triples
    #Disconnect graphs might become connected because of the nonterminal edges
    def allTriples(self, hgraph, var_mapping, ext_mapping, components=None, node_indexer=None):

        def checkChildren(components, index):
            ret = []
            remained = []
            for node in components:
                if index in node.frag.roots:
                    ret.append(node)
                else:
                    remained.append(node)
            return remained, ret

        stack = copy.copy(self.roots)
        stack = [(index, 0) for index in stack]
        stack.reverse()

        all_triples = []

        amr = self.graph

        visited_index = set()

        #attach_list = {}

        edge_alignment = bitarray(len(amr.edges))
        if edge_alignment.count() != 0:
            edge_alignment ^= edge_alignment

        if components:
            for node in components:
                if node.cut == 1:
                    edge_alignment |= node.frag.edges

        stack.reverse()
        while stack:
            curr_index, depth = stack.pop()
            if curr_index in visited_index:
                continue

            curr_node = amr.nodes[curr_index]
            if depth == 0:
                curr_ident = hgraph.registerNode(curr_index, var_mapping, ext_mapping)
                hgraph.roots.append(curr_ident)

            visited_index.add(curr_index)

            if components:
                components, nonterm_nodes = checkChildren(components, curr_index)
                for node in nonterm_nodes:
                    curr_frag = node.frag
                    suffix, other_indices = node.frag.extSuffix(curr_index)

                    curr_edge_label = 'A%d-%s' % (len(curr_frag.ext_list), suffix)
                    curr_edge_id = node_indexer[node]
                    nonterm_edge = NonterminalLabel(curr_edge_label, curr_edge_id)
                    all_triples.append((curr_index, nonterm_edge,
                        other_indices))

                    exts = [(index, depth+1) for index in other_indices]
                    stack.extend(exts)

            concept_edge = curr_node.c_edge
            if self.edges[concept_edge] == 1 and edge_alignment[concept_edge] != 1:
                child = None
                edge_label = str(curr_node)
                all_triples.append((curr_index, edge_label, child))

            for edge_index in curr_node.v_edges:
                curr_edge = amr.edges[edge_index]
                if self.edges[edge_index] == 1 and edge_alignment[edge_index] != 1:
                    tail_index = curr_edge.tail
                    all_triples.append((curr_index, curr_edge.label, tail_index))
                    stack.append((tail_index, depth+1))

        return all_triples

    def set_edges(self, edges):
        self.edges = edges

    def set_nodes(self, nodes):
        self.nodes = nodes

    def set_edge(self, edge_num):
        self.edges[edge_num] = 1

    def set_node(self, node_num):
        self.nodes[node_num] = 1

    def set_ext_set(self, ext_set):
        self.ext_set = ext_set

    def node_list(self):
        return [i for (i, existed) in enumerate(self.nodes) if existed]

    def edge_list(self):
        return [i for (i, existed) in enumerate(self.edges) if existed]

    def ext_nodes_str(self):
        s = ''
        for node_index in self.ext_set:
            s += str(self.graph.nodes[node_index])
            s += ' '
        return s.strip()

    def create_ext_mapping(self):
        num = 0
        for node_index in self.ext_set:
            self.ext_mapping[node_index] = '*%d' % num
            num += 1

    @staticmethod
    def initialize_from_alignment(nodes, edges, graph=None):
        frag = AMRFragment(len(edges), len(nodes), graph)
        frag.edges = edges
        frag.nodes = nodes
        return frag

    def __eq__(self, other):
        return ((self.edges == other.edges) and (self.nodes == other.nodes))

    def __hash__(self):
        s = ''
        if self.nodes != None:
            s += str(self.nodes)
        s += '-'
        if self.edges != None:
            s += str(self.edges)
        return hash(s)

    def is_ext(self, node_num):
        curr_node = self.graph.nodes[node_num]
        for edge_index in curr_node.edge_set():
            if self.edges[edge_index] == 0:
                return True
        return False

    def frag_repr(self, node_index, visited_index):
        assert node_index not in visited_index
        visited_index.add(node_index)

        result = '.'
        node = self.graph.nodes[node_index]
        const_edge_index = node.c_edge

        if self.edges[const_edge_index] == 1:
            result = '%s :%s' % (result, str(node))

        incoming_edges = []
        for edge_index in node.v_edges:
            if self.edges[edge_index] == 0:
                continue
            incoming_edges.append(edge_index)

        if len(incoming_edges) > 1:
            incoming_edges = sorted(incoming_edges, key=lambda index: str(self.graph.edges[index]))

        for edge_index in incoming_edges:
            curr_edge = self.graph.edges[edge_index]
            if curr_edge.tail in visited_index:
                tail_node = self.graph.nodes[curr_edge.tail]
                if self.edges[tail_node.c_edge] == 1:
                    result = '%s :%s (%s)' % (result, str(curr_edge), str(self.graph.edges[self.graph.nodes[curr_edge.tail].c_edge]))
                else:
                    result = '%s :%s .' % (result, str(curr_edge))
            else:
                sub_frag_repr = self.frag_repr(curr_edge.tail, visited_index)
                if ':' in sub_frag_repr:
                    sub_frag_repr = '(%s)' % sub_frag_repr
                result = '%s :%s %s' % (result, str(curr_edge), sub_frag_repr)
        return result

    def rootLabels(self):
        return [self.graph.nodes[node_index].node_str() for node_index in self.roots]

    def build_ext_set(self):
        self.ext_set = set(self.ext_list)

    #Perform dfs in a sorted way
    def build_ext_list(self):

        self.ext_list = []
        visited_index = set()

        stack = copy.copy(self.roots)
        stack.reverse()

        while stack:
            curr_node_index = stack.pop()
            assert self.nodes[curr_node_index] == 1
            if curr_node_index in visited_index: #Newly added, there might be reentrancy to same node between two nodes
                continue

            visited_index.add(curr_node_index)

            if self.is_ext(curr_node_index):
                self.ext_list.append(curr_node_index)

            curr_node = self.graph.nodes[curr_node_index]

            incoming_edges = []
            for edge_index in curr_node.v_edges:
                if self.edges[edge_index] == 0:
                    continue

                curr_edge = self.graph.edges[edge_index]

                incoming_edges.append((str(curr_edge), edge_index))

            if len(incoming_edges) > 1:
                incoming_edges = sorted(incoming_edges, key=lambda elem: elem[0])
            for (_, edge_index) in reversed(incoming_edges):
                curr_edge = self.graph.edges[edge_index]
                tail_index = curr_edge.tail
                stack.append(tail_index)

        for (i, ext_index) in enumerate(self.ext_list):
            self.ext_mapping[ext_index] = i

    # Compute the number of nodes connecting with this fragment.
    def build_outside_set(self, node_set, pi_map):
        for curr in self.ext_list:
            assert curr in node_set
            curr_node = self.graph.nodes[curr]
            for edge_index in curr_node.v_edges:
                tail_index = self.graph.edges[edge_index].tail
                assert self.graph.edges[edge_index].head == curr
                if tail_index not in node_set:
                    assert tail_index in pi_map
                    self.outside_set.add(pi_map[tail_index])
            for edge_index in curr_node.p_edges:
                head_index = self.graph.edges[edge_index].head
                assert self.graph.edges[edge_index].tail == curr
                if head_index not in node_set:
                    assert head_index in pi_map
                    self.outside_set.add(pi_map[head_index])

    def __str__(self):
        from util.hgraph.hgraph import Hgraph
        # s = Sample(None, 0)
        g = Hgraph()
        # frag_node = FragmentHGNode("DUMMY", self.start, self.end, self)
        # s.buildHgraph(g, frag_node, None, {}, {})
        return str(g)

    def missing_edges(self):
        s = ""
        for i in xrange(len(self.edges)):
            if self.edges[i] == 0:
                s += self.graph.edges[i].label
                s += '  '
        return s

    def var_from_graph(amr_graph, binary_reps):
        bits = binary_reps.split('.')
        assert bits[0] == '0', 'All binary representation should start from root 0'

#If two fragments are disjoint, then they do not share any common edge
def check_disjoint(f1, f2, capacities):
    result = f1.edges & f2.edges
    not_conflict = True
    for i in xrange(len(result)):
        if result[i] == 1 and capacities[i] <= 1:
            not_conflict = False
    return not_conflict

#Root operation is a bold guess: that the root of combination must be a root in one of the child fragment
def combine_fragments(f1, f2, refine=False):
    f1_rooted = (f2.root in f1.ext_set)
    f2_rooted = (f1.root in f2.ext_set)

    assert not refine
    #Currently we don't allow combination of subgraph having two roots
    if not (f1_rooted or f2_rooted):
        return None

    n_nodes = len(f1.nodes)
    n_edges = len(f1.edges)
    amr_graph = f1.graph
    new_frag = AMRFragment(n_edges, n_nodes, amr_graph)

    if f1.root != f2.root and f1_rooted and f2_rooted: #There is a cycle
        f1_rootnode = amr_graph.nodes[f1.root]
        f2_rootnode = amr_graph.nodes[f2.root]

        new_frag.root = f2.root if len(f2_rootnode.p_edges) > len(f1_rootnode.p_edges) else f1.root
    elif f1_rooted:
        new_frag.root = f1.root
    else:
        new_frag.root = f2.root
    nodes = f1.nodes | f2.nodes
    edges = f1.edges | f2.edges

    new_frag.set_edges(edges)
    new_frag.set_nodes(nodes)

    new_frag.build_ext_list()
    new_frag.build_ext_set()

    #Setting the word span of the new fragment
    if f1.start == -1:
        new_frag.setSpan(f2.start, f2.end)
    elif f2.start == -1:
        new_frag.setSpan(f1.start, f1.end)
    else:

        if f1.start < f2.start:
            assert f1.end <= f2.start, 'overlapping fragments'
            new_frag.set_span(f1.start, f2.end)
        else:
            new_frag.set_span(f2.start, f1.end)
    return new_frag

def find_unaligned_edge(curr_index, another_index, amr_graph, edge_alignment):
    curr_node = amr_graph.nodes[curr_index]

    n_nodes = len(amr_graph.nodes)
    n_edges = len(amr_graph.edges)
    unaligned_frag = AMRFragment(n_edges, n_nodes, amr_graph)

    for edge_index in curr_node.p_edges:
        if edge_alignment[edge_index] == 1:
            continue
        curr_pedge = amr_graph.edges[edge_index]
        p_index = curr_pedge.head

        if p_index == another_index:
            unaligned_frag.setRoot(p_index)

            unaligned_frag.set_node(p_index)
            unaligned_frag.set_node(curr_index)

            unaligned_frag.set_edge(edge_index)

            unaligned_frag.build_ext_list()
            unaligned_frag.build_ext_set()
            return unaligned_frag

    return None

def find_unaligned_path(curr_index, frag, edge_alignment, refine=False):
    amr_graph = frag.graph
    curr_node = amr_graph.nodes[curr_index]
    path = []

    n_nodes = len(amr_graph.nodes)
    n_edges = len(amr_graph.edges)
    unaligned_frag = AMRFragment(n_edges, n_nodes, amr_graph)

    for edge_index in curr_node.p_edges:
        if edge_alignment[edge_index] == 1:
            continue
        curr_pedge = amr_graph.edges[edge_index]
        p_index = curr_pedge.head

        #If two fragments are connected through one relation
        if p_index in frag.ext_set:
            ext_set = set()
            unaligned_frag.setRoot(p_index)

            unaligned_frag.set_node(p_index)
            unaligned_frag.set_node(curr_index)

            unaligned_frag.set_edge(edge_index)

            unaligned_frag.build_ext_list()
            unaligned_frag.build_ext_set()
            if refine:
                init_ext_frag(unaligned_frag)
            return unaligned_frag

    return None

#Combine two fragments, take union of the nodes
#Then traverse from each root and connect all the edges amongst the nodes
def general_combine_fragments(f1, f2):

    n_nodes = len(f1.nodes)
    n_edges = len(f1.edges)
    amr_graph = f1.graph
    new_frag = AMRFragment(n_edges, n_nodes, amr_graph)

    #connect_vertices = set()  #The nodes that connect the two fragments

    nodes = f1.nodes | f2.nodes
    edges = f1.edges | f2.edges

    new_frag.set_nodes(nodes)
    node_set = set(new_frag.node_list())
    # outside_set = set()
    for start_index in node_set:
        curr_node = amr_graph.nodes[start_index]
        for edge_index in curr_node.v_edges:
            curr_edge = amr_graph.edges[edge_index]
            if curr_edge.tail in node_set:
                if edges[edge_index] == 0:
                    edges[edge_index] = 1

    new_frag.set_edges(edges)

    roots = amr_graph.buildFragment(nodes)
    new_frag.roots = roots

    new_frag.build_ext_list()
    new_frag.build_ext_set()
    # new_frag.build_outside_set(node_set)

    if f1.start == -1:
        new_frag.setSpan(f2.start, f2.end)
    elif f2.start == -1:
        new_frag.setSpan(f1.start, f1.end)
    else:
        if f1.start < f2.start:
            assert f1.end <= f2.start, 'overlapping fragments'
            new_frag.setSpan(f1.start, f2.end)
        else:
            new_frag.setSpan(f2.start, f1.end)

    return new_frag, len(f1.ext_set | f2.ext_set)

def check_consist(parent, children):
    nodes = children[0].frag.nodes | children[1].frag.nodes
    edges = children[0].frag.edges | children[1].frag.edges
    for i in xrange(2, len(children)):
        nodes |= children[i].frag.nodes
        edges |= children[i].frag.edges
    if parent.frag.nodes != nodes or parent.frag.edges != edges:
        return False
    return True

#Initialize the external mapping for the initialize fragments
def init_ext_frag(frag, is_pred=False, is_op=False):

    #First initialize the root
    amr_graph = frag.graph
    root_node = amr_graph.nodes[frag.root]
    if frag.edges[root_node.c_edge] == 0:
        frag.ext_label[frag.root] = '.' #If a fragment does not have a root label, the root is identified as dot

    else:
        root_label = str(root_node)
        #if 'ARG' in str(frag) and ('of' not in str(frag)) and '-' in root_label: #Indicate we are dealing with a predicate structure
        if is_pred:
            frag.ext_label[frag.root] = 'PRED'
        elif is_op:
            frag.ext_label[frag.root] = root_label if '/' not in root_label else root_label.split('/')[1]
        else:
            frag.ext_label[frag.root] = 'ENT'

    for ext_index in frag.ext_set: #I assume only could be nodes without concept
        if ext_index == frag.root:
            continue

        curr_node = amr_graph.nodes[ext_index]
        assert frag.edges[curr_node.c_edge] == 0, 'weird node with concept found'

        #assert len(amr_graph.nodes[ext_index].p_edges) == 1, 'not only 1 parent found'
        for p_edge_index in curr_node.p_edges:
            if frag.edges[p_edge_index] == 1:
                curr_edge_label = str(amr_graph.edges[p_edge_index])
                frag.ext_label[ext_index] = curr_edge_label

def connect_adjacent(frags, logger):
    visited = set()
    new_frag_list = []
    updated = True
    for frag in frags:
        if frag in visited:
            continue
        visited.add(frag)
        curr_result = frag
        updated = False
        for cand_frag in frags:
            if cand_frag in visited:
                continue
            if check_disjoint(curr_result, cand_frag):

                new_result = combine_fragments(curr_result, cand_frag)
                if new_result:
                    updated = True
                    visited.add(cand_frag)
                    visited.add(curr_result)
                    curr_result = new_result
                else:
                    logger.writeln("Before:")
                    logger.writeln(str(cand_frag))
                    logger.writeln(str(curr_result))
                    logger.writeln("After:")
                    logger.writeln(str(new_result))

        new_frag_list.append(curr_result)
    return new_frag_list

