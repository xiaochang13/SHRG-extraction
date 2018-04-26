#!/usr/bin/python
import sys
#import time
import pickle
import os
import re
import cPickle
from amr_graph import *
from amr_utils import *
import logger
import argparse
from re_utils import *
from preprocess import *
from collections import defaultdict
from entities import identify_entities
from constants import *
from date_extraction import *
from cacheTransitionParser import *
import random
from sdp_process import *

def build_bimap(tok2frags):
    frag2map = defaultdict(set)
    index2frags = defaultdict(set)
    for index in tok2frags:
        for frag in tok2frags[index]:
            index2frags[index].add(frag)
            frag2map[frag].add(index)

    return (index2frags, frag2map)

#Here we try to make the tok to fragment mapping one to one
def rebuild_fragment_map(tok2frags):
    (index2frags, frag2map) = build_bimap(tok2frags)
    for index in tok2frags:
        if len(tok2frags[index]) > 1:
            new_frag_list = []
            min_frag = None
            min_length = 100
            for frag in tok2frags[index]:
                index_set = frag2map[frag]
                assert index in index_set
                if len(index_set) > 1:
                    if len(index_set) < min_length:
                        min_length = len(index_set)
                        min_frag = frag
                    index_set.remove(index)
                else:
                    new_frag_list.append(frag)
            if len(new_frag_list) == 0:
                assert min_frag is not None
                new_frag_list.append(min_frag)
            tok2frags[index] = new_frag_list
    return tok2frags

def mergeSpans(index_to_spans):
    new_index_to_spans = {}
    for index in index_to_spans:
        span_list = index_to_spans[index]
        span_list = sorted(span_list, key=lambda x:x[1])
        new_span_list = []
        curr_start = None
        curr_end = None

        for (idx, (start, end, _)) in enumerate(span_list):
            if curr_end is not None:
                #assert start >= curr_end, span_list
                if start < curr_end:
                    continue
                if start > curr_end: #There is a gap in between
                    new_span_list.append((curr_start, curr_end, None))
                    curr_start = start
                    curr_end = end
                else: #They equal, so update the end
                    curr_end = end

            else:
                curr_start = start
                curr_end = end

            if idx + 1 == len(span_list): #Have reached the last position
                new_span_list.append((curr_start, curr_end, None))
        new_index_to_spans[index] = new_span_list
    return new_index_to_spans

def getDateAttr(frag):
    date_relations = set(['time', 'year', 'month', 'day', 'weekday', 'century', 'era', 'decade', 'dayperiod', 'season', 'timezone'])
    root_index = frag.root
    amr_graph = frag.graph
    root_node = amr_graph.nodes[root_index]

    index_to_attr = {}

    attr_indices = set()
    for edge_index in root_node.v_edges:
        curr_edge = amr_graph.edges[edge_index]
        if curr_edge.label in date_relations:
            attr_indices.add(curr_edge.tail)
            index_to_attr[curr_edge.tail] = curr_edge.label.upper()
    return (attr_indices, index_to_attr)

#Given an alignment and the fragment, output the covered span
def getSpanSide(toks, alignments, frag, unaligned_toks):
    aligned_set = set()
    amr_graph = frag.graph

    covered_set = set()
    all_date_attrs, index_to_attr = getDateAttr(frag)

    index_to_toks = defaultdict(list)

    for curr_align in reversed(alignments):
        curr_tok = curr_align.split('-')[0]
        curr_frag = curr_align.split('-')[1]

        start = int(curr_tok)
        end = start + 1

        aligned_set.add(start)

        (index_type, index) = amr_graph.get_concept_relation(curr_frag)
        if index_type == 'c':
            if frag.nodes[index] == 1: #Covered current
                covered_set.add(start)
                index_to_toks[index].append(start)
                if index in all_date_attrs:
                    all_date_attrs.remove(index)

        else: #An edge covered span
            if frag.edges[index] == 1:
                covered_set.add(start)

    covered_toks = sorted(list(covered_set))
    non_covered = [amr_graph.nodes[index].node_str() for index in all_date_attrs]
    return covered_toks, non_covered, index_to_toks

def extractNodeMapping(alignments, amr_graph):
    aligned_set = set()

    node_to_span = defaultdict(list)
    edge_to_span = defaultdict(list)

    num_nodes = len(amr_graph.nodes)
    num_edges = len(amr_graph.edges)

    op_toks = []
    role_toks = []
    for curr_align in reversed(alignments):
        curr_tok = curr_align.split('-')[0]
        curr_frag = curr_align.split('-')[1]

        start = int(curr_tok)
        end = start + 1

        aligned_set.add(start)

        (index_type, index) = amr_graph.get_concept_relation(curr_frag)
        if index_type == None:
            return (None, None, None, None, None)

        if index_type == 'c':
            node_to_span[index].append((start, end, None))
            curr_node = amr_graph.nodes[index]

            #Extract ops for entities
            if len(curr_node.p_edges) == 1:
                par_edge = amr_graph.edges[curr_node.p_edges[0]]
                if 'op' == par_edge.label[:2]:
                    op_toks.append((start, curr_node.c_edge))

            if curr_node.is_named_entity():
                role_toks.append((start, curr_node.c_edge))

        else:
            edge_to_span[index].append((start, end, None))
    new_node_to_span = mergeSpans(node_to_span)
    new_edge_to_span = mergeSpans(edge_to_span)

    return (op_toks, role_toks, new_node_to_span, new_edge_to_span, aligned_set)

def extract_fragments(alignments, amr_graph):
    tok2frags = defaultdict(list)

    num_nodes = len(amr_graph.nodes)
    num_edges = len(amr_graph.edges)

    op_toks = []
    role_toks = []
    for curr_align in reversed(alignments):
        curr_tok = curr_align.split('-')[0]
        curr_frag = curr_align.split('-')[1]

        start = int(curr_tok)
        end = start + 1

        (index_type, index) = amr_graph.get_concept_relation(curr_frag)
        frag = AMRFragment(num_edges, num_nodes, amr_graph)
        if index_type == 'c':
            frag.set_root(index)
            curr_node = amr_graph.nodes[index]

            #Extract ops for entities
            if len(curr_node.p_edges) == 1:
                par_edge = amr_graph.edges[curr_node.p_edges[0]]
                if 'op' == par_edge.label[:2]:
                    op_toks.append((start, curr_node.c_edge))

            if curr_node.is_named_entity():
                role_toks.append((start, curr_node.c_edge))

            frag.set_edge(curr_node.c_edge)

        else:
            frag.set_edge(index)
            curr_edge = amr_graph.edges[index]
            frag.set_root(curr_edge.head)
            frag.set_node(curr_edge.tail)

        frag.build_ext_list()
        frag.build_ext_set()

        tok2frags[start].append(frag)

    for index in tok2frags:
        if len(tok2frags[index]) > 1:
            tok2frags[index] = connect_adjacent(tok2frags[index], logger)

    tok2frags = rebuild_fragment_map(tok2frags)
    for index in tok2frags:
        for frag in tok2frags[index]:
            frag.set_span(index, index+1)

    return (op_toks, role_toks, tok2frags)

#Verify this fragment contains only one edge and return it
def unique_edge(frag):
    #assert frag.edges.count() == 1, 'Not unify edge fragment found'
    amr_graph = frag.graph
    edge_list = []
    n_edges = len(frag.edges)
    for i in xrange(n_edges):
        if frag.edges[i] == 1:
            edge_list.append(i)
    assert len(edge_list) == frag.edges.count()
    return tuple(edge_list)

class AMR_stats(object):
    def __init__(self):
        self.num_reentrancy = 0
        self.num_predicates = defaultdict(int)
        self.num_nonpredicate_vals = defaultdict(int)
        self.num_consts = defaultdict(int)
        self.num_named_entities = defaultdict(int)
        self.num_entities = defaultdict(int)
        self.num_relations = defaultdict(int)

    def update(self, local_re, local_pre, local_non, local_con, local_ent, local_ne):
        self.num_reentrancy += local_re
        for s in local_pre:
            self.num_predicates[s] += local_pre[s]

        for s in local_non:
            self.num_nonpredicate_vals[s] += local_non[s]

        for s in local_con:
            self.num_consts[s] += local_con[s]

        for s in local_ent:
            self.num_entities[s] += local_ent[s]

        for s in local_ne:
            self.num_named_entities[s] += local_ne[s]
        #for s in local_rel:
        #    self.num_relations[s] += local_rel[s]

    def collect_stats(self, amr_graphs):
        for amr in amr_graphs:
            (named_entity_nums, entity_nums, predicate_nums, variable_nums, const_nums, reentrancy_nums) = amr.statistics()
            self.update(reentrancy_nums, predicate_nums, variable_nums, const_nums, entity_nums, named_entity_nums)

    def dump2dir(self, dir):
        def dump_file(f, dict):
            sorted_dict = sorted(dict.items(), key=lambda k:(-k[1], k[0]))
            for (item, count) in sorted_dict:
                print >>f, '%s %d' % (item, count)
            f.close()

        pred_f = open(os.path.join(dir, 'pred'), 'w')
        non_pred_f = open(os.path.join(dir, 'non_pred_val'), 'w')
        const_f = open(os.path.join(dir, 'const'), 'w')
        entity_f = open(os.path.join(dir, 'entities'), 'w')
        named_entity_f = open(os.path.join(dir, 'named_entities'), 'w')
        #relation_f = open(os.path.join(dir, 'relations'), 'w')

        dump_file(pred_f, self.num_predicates)
        dump_file(non_pred_f, self.num_nonpredicate_vals)
        dump_file(const_f, self.num_consts)
        dump_file(entity_f, self.num_entities)
        dump_file(named_entity_f, self.num_named_entities)
        #dump_file(relation_f, self.num_relations)

    def loadFromDir(self, dir):
        def load_file(f, dict):
            for line in f:
                item = line.strip().split(' ')[0]
                count = int(line.strip().split(' ')[1])
                dict[item] = count
            f.close()

        pred_f = open(os.path.join(dir, 'pred'), 'r')
        non_pred_f = open(os.path.join(dir, 'non_pred_val'), 'r')
        const_f = open(os.path.join(dir, 'const'), 'r')
        entity_f = open(os.path.join(dir, 'entities'), 'r')
        named_entity_f = open(os.path.join(dir, 'named_entities'), 'r')

        load_file(pred_f, self.num_predicates)
        load_file(non_pred_f, self.num_nonpredicate_vals)
        load_file(const_f, self.num_consts)
        load_file(entity_f, self.num_entities)
        load_file(named_entity_f, self.num_named_entities)

    def __str__(self):
        s = ''
        s += 'Total number of reentrancies: %d\n' % self.num_reentrancy
        s += 'Total number of predicates: %d\n' % len(self.num_predicates)
        s += 'Total number of non predicates variables: %d\n' % len(self.num_nonpredicate_vals)
        s += 'Total number of constants: %d\n' % len(self.num_consts)
        s += 'Total number of entities: %d\n' % len(self.num_entities)
        s += 'Total number of named entities: %d\n' % len(self.num_named_entities)

        return s

def similarTok(curr_var, tok):
    if curr_var == tok:
        return True
    var_len = len(curr_var)
    tok_len = len(tok)
    if var_len > 3 and tok_len > 3 and tok[:4] == curr_var[:4]:
        return True
    if isNum(tok) and tok in curr_var:
        return True
    if isNum(curr_var) and curr_var in tok:
        return True

#Traverse AMR from top down, also categorize the sequence in case of alignment existed
def depthFirstPi(amr, tok_seq, all_alignments, sorted_indexes, print_info=False):

    old_depth = -1
    depth = -1
    stack = [(amr.root, TOP, None, False)] #Start from the root of the AMR
    aux_stack = []
    seq = []

    cate_tok_seq = []
    ret_index = 0

    seq_map = {}   #Map each span to a category
    visited = set()

    covered = set()
    multiple_covered = set()

    node_to_label = {}
    cate_to_index = {}

    node_to_parents = defaultdict(list)  #A map from each aligned index to its unaligned parent(s)
    adjacent_nodes = defaultdict(set) #Build a map from node to all its adjacent nodes

    sorted_set = set(sorted_indexes)
    verb_map = {}
    pred_freq_thre = 0
    var_freq_thre = 0

    edge_list = []
    edge_map = defaultdict(set)
    topo_list = []

    while stack:
        curr_node_index, rel, parent_index, unaligned = stack.pop()

        curr_node = amr.nodes[curr_node_index]
        curr_var = curr_node.node_str()

        i = 0

        if curr_node_index in visited: #A reentrancy found
            continue

        topo_list.append(curr_node_index)
        visited.add(curr_node_index)

        curr_node = amr.nodes[curr_node_index]
        if amr.is_named_entity(curr_node):
            exclude_rels = ['wiki', 'name']
        elif amr.is_date_entity(curr_node):
            exclude_rels = date_relations
        else:
            exclude_rels = []

        #unaligned = curr_node_index not in sorted_set
        for edge_index in reversed(curr_node.v_edges):
            curr_edge = amr.edges[edge_index]
            child_index = curr_edge.tail
            if curr_edge.label in exclude_rels:
                visited.add(child_index)
                if curr_edge.label != 'name':
                    tail_node = amr.nodes[child_index]
                    for next_edge_index in reversed(tail_node.v_edges):
                        next_edge = amr.edges[next_edge_index]
                        next_child_index = next_edge.tail
                        edge_map[curr_node_index].add(next_child_index)
                        edge_map[next_child_index].add(curr_node_index)
                        edge_list.append((curr_node_index, child_index))

                        stack.append((next_child_index, next_edge.label, curr_node_index, unaligned))
                continue

            edge_map[curr_node_index].add(child_index)
            edge_map[child_index].add(curr_node_index)
            edge_list.append((curr_node_index, child_index))

            stack.append((child_index, curr_edge.label, curr_node_index, unaligned))

    topo_set = set(topo_list)

    assert len(topo_set) == len(topo_list)
    for node_index in topo_list:
        #assert node_index in edge_map
        if node_index not in edge_map:
            print str(amr)
        #if node_index not in edge_map:
        #    print 'Ignored edge here:'
        #    print amr.nodes[node_index].node_str()
        #    continue

    if print_info:
        print topo_list
        print edge_map
        print edge_list

    return topo_list, edge_map, edge_list


#Traverse AMR from top down, also categorize the sequence in case of alignment existed
def buildPiSeq(amr, tok_seq, all_alignments, sorted_indexes, print_info=False):

    old_depth = -1
    depth = -1
    stack = [(amr.root, TOP, None, False)] #Start from the root of the AMR
    aux_stack = []
    seq = []

    cate_tok_seq = []
    ret_index = 0

    seq_map = {}   #Map each span to a category
    visited = set()

    covered = set()
    multiple_covered = set()

    node_to_label = {}
    cate_to_index = {}

    node_to_parents = defaultdict(list)  #A map from each aligned index to its unaligned parent(s)
    adjacent_nodes = defaultdict(set) #Build a map from node to all its adjacent nodes

    sorted_set = set(sorted_indexes)
    verb_map = {}
    pred_freq_thre = 0
    var_freq_thre = 0

    vertex_set = set()

    edge_list = []
    edge_map = defaultdict(set)

    while stack:
        curr_node_index, rel, parent_index, unaligned = stack.pop()
        vertex_set.add(curr_node_index)
        if unaligned:
            node_to_parents[curr_node_index].append(parent_index)

        curr_node = amr.nodes[curr_node_index]
        curr_var = curr_node.node_str()

        i = 0

        if curr_node_index in visited: #A reentrancy found
            continue

        visited.add(curr_node_index)

        if curr_node_index in all_alignments:
            try:
                assert curr_node_index in sorted_set
            except:
                print 'weird here'
                print sorted_indexes
                print all_alignments.keys()
                print amr.nodes[curr_node_index].node_str()
                print str(amr)
            exclude_rels, cur_symbol, categorized = amr.get_symbol(curr_node_index, verb_map, pred_freq_thre, var_freq_thre)

        else:
            exclude_rels, cur_symbol, categorized = [], curr_var, False

        unaligned = curr_node_index not in sorted_set
        for edge_index in reversed(curr_node.v_edges):
            curr_edge = amr.edges[edge_index]
            child_index = curr_edge.tail
            if curr_edge.label in exclude_rels:
                visited.add(child_index)
                #if 'VERBAL' in cur_symbol: #Might have other relations
                if curr_edge.label != 'name':
                    #assert False
                    tail_node = amr.nodes[child_index]
                    for next_edge_index in reversed(tail_node.v_edges):
                        next_edge = amr.edges[next_edge_index]
                        next_child_index = next_edge.tail
                        edge_map[curr_node_index].add(next_child_index)
                        edge_map[next_child_index].add(curr_node_index)
                        edge_list.append((curr_node_index, child_index))

                        stack.append((next_child_index, next_edge.label, curr_node_index, unaligned))
                continue

            edge_map[curr_node_index].add(child_index)
            edge_map[child_index].add(curr_node_index)
            edge_list.append((curr_node_index, child_index))

            stack.append((child_index, curr_edge.label, curr_node_index, unaligned))

    traversed = set()
    topo_list = []
    appeared = set()
    for node_index in sorted_indexes:
        if node_index in traversed:
            continue

        if node_index not in edge_map:
            print 'Ignored edge here:'
            print amr.nodes[node_index].node_str()
            continue

        traversed.add(node_index)
        unaligned_seq = topSortUnaligned(node_index, node_to_parents, traversed)
        for index in unaligned_seq:
            if index not in appeared:
                topo_list.append(index)
        topo_list.append(node_index)

    unvisited = vertex_set - set(topo_list)
    for vertex in unvisited:
        topo_list.append(vertex)

    if print_info:
        print sorted_indexes
        print all_alignments.keys()
        print vertex_set
        print topo_list

    if len(vertex_set) != len(topo_list):
        uncovered = set(topo_list) - vertex_set
        for index in uncovered:
            print amr.nodes[index].node_str()
    #assert len(vertex_set) == len(topo_list)
    return topo_list, edge_map, edge_list

#def deterministicTreeDecomposition(edge_triples, pi_seq):
def topSortUnaligned(node_index, node_to_parents, visited):
    topo_list = []
    stack = [node_index]
    while stack:
        curr_index = stack.pop()
        if curr_index not in visited:
            topo_list.append(curr_index)
        visited.add(curr_index)
        if curr_index in node_to_parents:
            for p_index in node_to_parents[curr_index]:
                if p_index not in visited:
                    stack.append(p_index)
    return topo_list[::-1]

def categorizedSpans(all_alignments, node_to_label):
    visited = set()

    all_alignments = sorted(all_alignments.items(), key=lambda x:len(x[1]))
    span_map = {}
    end_index_map = {}

    for (node_index, aligned_spans) in all_alignments:
        if node_index in node_to_label:
            aligned_label = node_to_label[node_index]
            for (start, end, wiki_label) in aligned_spans:
                span_set = set(xrange(start, end))
                if len(span_set & visited) != 0:
                    continue

                visited |= span_set

                span_map[(start, end)] = (node_index, aligned_label, wiki_label)
                end_index_map[start] = end

    return span_map, end_index_map, visited

def orderedTreewidth(args):
    if not os.path.isdir(args.run_dir):
        os.system('mkdir -p %s' % args.run_dir)
    logger.file = open(os.path.join(args.run_dir, 'logger'), 'w')

    amr_file = os.path.join(args.data_dir, 'amr')
    alignment_file = os.path.join(args.data_dir, 'alignment')
    if args.use_lemma:
        tok_file = os.path.join(args.data_dir, 'lemmatized_token')
    else:
        tok_file = os.path.join(args.data_dir, 'token')
    pos_file = os.path.join(args.data_dir, 'pos')

    amr_graphs = load_amr_graphs(amr_file)
    alignments = [line.strip().split() for line in open(alignment_file, 'r')]
    toks = [line.strip().split() for line in open(tok_file, 'r')]

    poss = [line.strip().split() for line in open(pos_file, 'r')]

    assert len(amr_graphs) == len(alignments) and len(amr_graphs) == len(toks) and len(amr_graphs) == len(poss), '%d %d %d %d %d' % (len(amr_graphs), len(alignments), len(toks), len(poss))

    num_self_cycle = 0
    used_sents = 0

    amr_statistics = AMR_stats()

    if args.use_stats:
        amr_statistics.loadFromDir(args.stats_dir)
    else:
        os.system('mkdir -p %s' % args.stats_dir)
        amr_statistics.collect_stats(amr_graphs)
        amr_statistics.dump2dir(args.stats_dir)

    singleton_num = 0.0
    multiple_num = 0.0
    total_num = 0.0
    empty_num = 0.0

    aligned_indexes = set()

    width_dist = defaultdict(int)
    total_width = 0.0
    total_sents = 0.0

    all_piseqs = []
    all_edgemaps = []

    for (sent_index, (tok_seq, pos_seq, alignment_seq, amr)) in enumerate(zip(toks, poss, alignments, amr_graphs)):

        logger.writeln('Sentence #%d' % (sent_index+1))
        logger.writeln(' '.join(tok_seq))

        amr.setStats(amr_statistics)

        edge_alignment = bitarray(len(amr.edges))
        if edge_alignment.count() != 0:
            edge_alignment ^= edge_alignment
        assert edge_alignment.count() == 0

        concept_alignment = bitarray(len(amr.nodes))
        if concept_alignment.count() != 0:
            concept_alignment ^= concept_alignment
        assert concept_alignment.count() == 0

        has_cycle = False
        if amr.check_self_cycle():
            num_self_cycle += 1
            has_cycle = True

        amr.set_sentence(tok_seq)
        amr.set_poss(pos_seq)

        aligned_fragments = []
        reentrancies = {}  #Map multiple spans as reentrancies, keeping only one as original, others as connections

        has_multiple = False
        no_alignment = False

        aligned_set = set()

        (opt_toks, role_toks, node_to_span, edge_to_span, temp_aligned) = extractNodeMapping(alignment_seq, amr)
        if edge_to_span is None:
            print sent_index
            print str(amr)
            print ' '.join(tok_seq)
            print ' '.join(alignment_seq)
            continue

        temp_unaligned = set(xrange(len(pos_seq))) - temp_aligned

        all_frags = []
        all_alignments = defaultdict(list)

        ####Extract named entities#####
        for (frag, wiki_label) in amr.extract_entities():
            if len(opt_toks) == 0:
                logger.writeln("No alignment for the entity found")

            (aligned_indexes, entity_spans) = all_aligned_spans(frag, opt_toks, role_toks, temp_unaligned)
            root_node = amr.nodes[frag.root]

            entity_mention_toks = root_node.namedEntityMention()

            total_num += 1.0
            if entity_spans:
                concept_alignment |= frag.nodes
                entity_spans = removeRedundant(tok_seq, entity_spans, entity_mention_toks)
                if len(entity_spans) == 1:
                    singleton_num += 1.0
                    logger.writeln('Single fragment')
                    for (start, end) in entity_spans:
                        logger.writeln(' '.join(tok_seq[start:end]))
                        all_alignments[frag.root].append((start, end, wiki_label))
                        temp_aligned |= set(xrange(start, end))
                else:
                    multiple_num += 1.0
                    logger.writeln('Multiple fragment')
                    logger.writeln(aligned_indexes)
                    logger.writeln(' '.join([tok_seq[index] for index in aligned_indexes]))

                    for (start, end) in entity_spans:
                        logger.writeln(' '.join(tok_seq[start:end]))
                        all_alignments[frag.root].append((start, end, wiki_label))
                        temp_aligned |= set(xrange(start, end))
            else:
                empty_num += 1.0

        ####Process date entities
        date_entity_frags = amr.extract_all_dates()
        for frag in date_entity_frags:
            all_date_indices, index_to_attr = getDateAttr(frag)
            covered_toks, non_covered, index_to_toks = getSpanSide(tok_seq, alignment_seq, frag, temp_unaligned)

            covered_set = set(covered_toks)

            all_spans = getContinuousSpans(covered_toks, temp_unaligned, covered_set)
            if all_spans:
                concept_alignment |= frag.nodes
                temp_spans = []
                for start, end in all_spans:
                    if start > 0 and (start-1) in temp_unaligned:
                        if tok_seq[start-1] in str(frag) and tok_seq[0] in '0123456789':
                            temp_spans.append((start-1, end))
                        else:
                            temp_spans.append((start, end))
                    else:
                        temp_spans.append((start, end))
                all_spans = temp_spans
                all_spans = removeDateRedundant(all_spans)
                for start, end in all_spans:
                    all_alignments[frag.root].append((start, end, None))
                    temp_aligned |= set(xrange(start, end))

            else:
                for index in temp_unaligned:
                    curr_tok = tok_seq[index]
                    found = False
                    for un_tok in non_covered:
                        if curr_tok[0] in '0123456789' and curr_tok in un_tok:
                            found = True
                            break
                    if found:
                        all_alignments[frag.root].append((index, index+1, None))
                        concept_alignment |= frag.nodes
                        temp_aligned.add(index)

        #Verbalization list
        #verb_map = {}
        #for (index, curr_tok) in enumerate(tok_seq):
        #    if curr_tok in VERB_LIST:

        #        for subgraph in VERB_LIST[curr_tok]:

        #            matched_frags = amr.matchSubgraph(subgraph)
        #            if matched_frags:
        #                temp_aligned.add(index)

        #            for (node_index, ex_rels) in matched_frags:
        #                all_alignments[node_index].append((index, index+1, None))
        #                verb_map[node_index] = subgraph

        for node_index in node_to_span:
            if node_index in all_alignments:
                continue

            if concept_alignment[node_index] == 1:
                continue

            if len(node_to_span[node_index]) > 0:
                all_alignments[node_index] = node_to_span[node_index]

        span_indexes = []
        for node_index in all_alignments:
            for (start, end, _) in all_alignments[node_index]:
                span_indexes.append((start, end, node_index))
        sorted_indexes = sorted(span_indexes, key=lambda x: (x[0], x[1]))
        sent_alignments = sorted_indexes
        #for (x, y, z) in sorted_indexes:
        #    print x, y, amr.nodes[z].node_str()
        sorted_indexes = [z for (x, y, z) in sorted_indexes]

        print_info = False
        if args.depth:
            pi_seq, edge_map, edge_list = depthFirstPi(amr, tok_seq, all_alignments, sorted_indexes, print_info)
        else:
            pi_seq, edge_map, edge_list = buildPiSeq(amr, tok_seq, all_alignments, sorted_indexes, print_info)

            if args.random:
                random.shuffle(pi_seq)
            elif args.reversed:
                pi_seq = pi_seq[::-1]

        all_piseqs.append(pi_seq)
        all_edgemaps.append(edge_map)

        #Next compute the distance table of the connection
        #edge_map_list = defaultdict(deque)

        #error_occur = False

    runMachine(toks, all_piseqs, all_edgemaps, amr_graphs)

def instanceTreewidth(pi_seq, edge_map):
    edge_map_list = defaultdict(deque)

    for (i, curr_vertex) in enumerate(pi_seq):
        for (j, next_vertex) in enumerate(pi_seq):
            if j <= i:
                continue
            if next_vertex in edge_map[curr_vertex]:
                curr_dist = len(pi_seq) - j
                edge_map_list[curr_vertex].append(curr_dist)
    dm = DeterministicMachine()
    curr_width = dm.start(pi_seq, edge_map_list, edge_map)
    return curr_width

def runMachine(all_toks, all_piseqs, all_edgemaps):
    n_instance = len(all_piseqs)
    width_dist = defaultdict(int)
    total_width = 0.0
    total_sents = 0.0
    for sent_index in xrange(n_instance):

        pi_seq = all_piseqs[sent_index]
        edge_map = all_edgemaps[sent_index]

        curr_width = instanceTreewidth(pi_seq, edge_map)

        total_width += float(curr_width)
        total_sents += 1.0
        width_dist[curr_width] += 1

    sorted_width = sorted(width_dist.items(), key=lambda x: -x[1])
    for (width, count) in sorted_width:
        print 'width: %d, count: %d' % (width, count)
    print 'Average width', total_width/total_sents
    #high_width_f.close()

def buildLinearEnt(entity_name, ops):
    ops_strs = ['op%d( %s )op%d' % (index, s, index) for (index, s) in enumerate(ops, 1)]
    ent_repr = '%s name( name %s )name' % (entity_name, ' '.join(ops_strs))
    return ent_repr

def isSpecial(symbol):
    for l in ['ENT', 'NE', 'VERB', 'SURF', 'CONST', 'DATE', 'VERBAL']:
        if l in symbol:
            return True
    return False

def getIndexedForm(linearized_tokseq):
    new_seq = []
    indexer = {}
    map_repr = []
    #cate_to_node_repr = {}
    for (start, end, (tok, node_repr, wiki_label)) in linearized_tokseq:
        if isSpecial(tok):
            new_tok = '%s-%d' % (tok, indexer.setdefault(tok, 0))
            indexer[tok] += 1
            new_seq.append(new_tok)
            #cate_to_node_repr[new_tok] = (node_repr, wiki_label)
            map_repr.append('%s++%d++%d++%s++%s' % (new_tok, start, end, node_repr, wiki_label))
        else:
            new_seq.append(tok)
    return new_seq, map_repr

#Given dev or test data, build the linearized token sequence
#Based on entity mapping from training, NER tagger
def conceptID(args):
    return

def dateMap(dateFile):
    dates_in_lines = []
    for line in open(dateFile):
        date_spans = []
        if line.strip():
            spans = line.strip().split()
            for sp in spans:
                start = int(sp.split('-')[0])
                end = int(sp.split('-')[1])
                date_spans.append((start, end))
        dates_in_lines.append(date_spans)
    return dates_in_lines

#Build the entity map for concept identification
#Choose either the most probable category or the most probable node repr
def loadMap(map_file):
    span_to_cate = {}

    #First load all possible mappings each span has
    with open(map_file, 'r') as map_f:
        for line in map_f:
            if line.strip():
                spans = line.strip().split('##')
                for s in spans:
                    try:
                        fields = s.split('++')
                        toks = fields[1]
                        wiki_label = fields[2]
                        node_repr = fields[3]
                        category = fields[-1]
                    except:
                        print spans, line
                        print fields
                        sys.exit(1)
                    if toks not in span_to_cate:
                        span_to_cate[toks] = defaultdict(int)
                    span_to_cate[toks][(category, node_repr, wiki_label)] += 1

    mle_map = {}
    for toks in span_to_cate:
        sorted_types = sorted(span_to_cate[toks].items(), key=lambda x:-x[1])
        mle_map[toks] = sorted_types[0][0]
    return mle_map

def semanticWidth(args):
    all_toks, all_piseqs, all_edgemaps = processData(args.input_file)
    runMachine(all_toks, all_piseqs, all_edgemaps)

#For each sentence, rebuild the map from categorized form to graph side nodes
def rebuildMap(args):
    return

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()

    argparser.add_argument("--input_file", type=str, help="the original AMR graph files", required=False)
    argparser.add_argument("--high_width_file", type=str, help="the original AMR graph files", required=False)
    argparser.add_argument("--stop", type=str, help="stop words file", required=False)
    argparser.add_argument("--lemma", type=str, help="lemma file", required=False)
    argparser.add_argument("--data_dir", type=str, help="the data directory for dumped AMR graph objects, alignment and tokenized sentences")
    argparser.add_argument("--map_file", type=str, help="map file from training")
    argparser.add_argument("--date_file", type=str, help="all the date in each sentence")
    argparser.add_argument("--run_dir", type=str, help="the output directory for saving the constructed forest")
    argparser.add_argument("--use_lemma", action="store_true", help="if use lemmatized tokens")
    argparser.add_argument("--random", action="store_true", help="if use random order")
    argparser.add_argument("--reversed", action="store_true", help="if to reverse the pi sequence")
    argparser.add_argument("--depth", action="store_true", help="if to reverse the pi sequence")
    argparser.add_argument("--parallel", action="store_true", help="if to linearize parallel sequences")
    argparser.add_argument("--use_stats", action="store_true", help="if use a built-up statistics")
    argparser.add_argument("--stats_dir", type=str, help="the statistics directory")
    argparser.add_argument("--min_prd_freq", type=int, default=50, help="threshold for filtering predicates")
    argparser.add_argument("--min_var_freq", type=int, default=50, help="threshold for filtering non predicate variables")
    argparser.add_argument("--index_unknown", action="store_true", help="if to index the unknown predicates or non predicate variables")

    args = argparser.parse_args()
    #semanticWidth(args)
    orderedTreewidth(args)
