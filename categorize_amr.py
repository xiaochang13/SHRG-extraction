#!/usr/bin/python
from amr_utils import *
from HRGSample import *
import argparse
from re_utils import *
from filter_stop_words import *
from preprocess import *
from date_extraction import *
from treewidth import depthFirstPi, buildPiSeq
import alignment_utils

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
    assert len(frag.roots) == 1
    root_index = frag.roots[0]
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

def collapseTokens(tok_seq, lemma_seq, pos_seq, span_to_type, isTrain=True):
    n_toks = len(tok_seq)
    collapsed = set()

    new_alignment = defaultdict(list)
    collapsed_seq = []
    collapsed_lem = []
    collapsed_pos = []
    for i in xrange(n_toks):
        if i in collapsed:
            continue
        for j in xrange(i+1, n_toks+1):
            if (i, j) in span_to_type:
                node_idx, subgraph_repr, curr_sym = span_to_type[(i, j)]
                collapsed |= set(xrange(i, j))
                curr_idx = len(collapsed_seq)
                new_alignment[node_idx].append((curr_idx, curr_idx+1, None))
                if 'NE_' in curr_sym or 'DATE' in curr_sym or 'NUMBER' in curr_sym or curr_sym == "NE":
                    if "NE_" in curr_sym:
                        rd = random.random()
                        if rd >= 0.9 and isTrain:
                            curr_sym = "NE"
                    collapsed_seq.append(curr_sym)
                    collapsed_lem.append(curr_sym)
                    if 'NE' in curr_sym:
                        collapsed_pos.append('NE')
                    elif 'DATE' in curr_sym:
                        collapsed_pos.append('DATE')
                    else:
                        collapsed_pos.append('NUMBER')
                elif j - i > 1:
                    collapsed_seq.append('@'.join(tok_seq[i:j]).lower())
                    collapsed_lem.append('@'.join(lemma_seq[i:j]).lower())
                    collapsed_pos.append('PHRASE')
                else:

                    collapsed_seq.append(tok_seq[j-1].lower())
                    collapsed_lem.append(lemma_seq[j-1].lower())
                    collapsed_pos.append(pos_seq[j-1])

        if i not in collapsed:
            if "LRB" in tok_seq[i]:
                tok_seq[i] = '('
                lemma_seq[i] = '('
            elif "RRB" in tok_seq[i]:
                tok_seq[i] = ')'
                lemma_seq[i] = ')'
            collapsed_seq.append(tok_seq[i].lower())
            collapsed_lem.append(lemma_seq[i].lower())
            collapsed_pos.append(pos_seq[i])

    return collapsed_seq, collapsed_lem, collapsed_pos, new_alignment

#Extract the span to fragment alignment, don't allow one to multiple, but one to multiple
def initAlignments(amr, tok_seq, pos_seq, all_alignments, unaligned, verb_map, nodeid_to_frag,
                   pred_freq_thre=50, var_freq_thre=50, use_pos=False):

    span_to_frag = {}
    depth = -1
    stack = [(amr.root, TOP, None, 0)] #Start from the root of the AMR
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

    while stack:
        old_depth = depth
        curr_node_index, incoming_edge_index, parent, depth = stack.pop()
        curr_node = amr.nodes[curr_node_index]
        curr_var = curr_node.node_str()

        curr_frag = None
        rel = amr.edges[incoming_edge_index].label if isinstance(incoming_edge_index, int) \
            else incoming_edge_index

        if curr_node_index in visited: #A reentrancy found
            seq.append((rel+LBR, None))
            seq.append((RET + ('-%d' % ret_index), None))
            ret_index += 1
            aux_stack.append((RBR+rel, None))
            continue

        visited.add(curr_node_index)
        seq.append((rel+LBR, None))

        if curr_node_index in all_alignments:
            curr_frag, exclude_rels, cur_symbol, categorized = amr.getFragment(
                curr_node_index, verb_map, pred_freq_thre, var_freq_thre, nodeid_to_frag)
            spans = all_alignments[curr_node_index]
            assert len(spans) > 0
            try:
                for start, end in spans:
                    span_to_frag[(start, end)] = curr_frag
            except:
                print all_alignments
                print span_to_frag
                print spans
                sys.exit(1)
        else:
            exclude_rels = []

        for edge_index in reversed(curr_node.v_edges):
            curr_edge = amr.edges[edge_index]
            child_index = curr_edge.tail
            if curr_edge.label in exclude_rels:
                visited.add(child_index)
                if cur_symbol[:2] != 'NE': #Might have other relations
                    tail_node = amr.nodes[child_index]
                    for next_edge_index in reversed(tail_node.v_edges):
                        next_edge = amr.edges[next_edge_index]
                        next_child_index = next_edge.tail
                        stack.append((next_child_index, next_edge.label, curr_var, depth+1))
                continue
            stack.append((child_index, curr_edge.label, curr_var, depth+1))

    return span_to_frag

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

def initialize_edge_alignment(aligned_fragments, edge_alignment):
    for frag in aligned_fragments:
        edge_alignment |= frag.edges

#This method output all the unaligned node information of the current AMR graph
def output_all_unaligned_nodes(edge_alignment, amr_graph):
    un_seq = []
    un_nodes = []
    for i in xrange(len(amr_graph.nodes)):
        curr_node = amr_graph.nodes[i]
        c_edge_index = curr_node.c_edge
        if edge_alignment[c_edge_index] == 0: #Found a concept that is not aligned
            un_seq.append(str(curr_node))
            un_nodes.append(curr_node)
    #print >> unalign_f, ' '.join(un_seq)
    return un_nodes

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

def build_bimap(tok2frags):
    frag2map = defaultdict(set)
    index2frags = defaultdict(set)
    for index in tok2frags:
        for frag in tok2frags[index]:
            index2frags[index].add(frag)
            frag2map[frag].add(index)
            #matched_list = extract_patterns(str(frag), '~e\.[0-9]+(,[0-9]+)*')
            #matched_indexes = parse_indexes(matched_list)
            #for matched_index in matched_indexes:
            #    frag2map[frag].add(matched_index)
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

def extract_fragments(s2g_alignment, amr_graph):
    alignments = s2g_alignment.strip().split()
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
                if 'op' in par_edge.label:
                    op_toks.append((start, curr_node.c_edge))

            if curr_node.is_entity():
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
            frag.setSpan(index, index+1)

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

def removeAligned(spans, aligned_toks):
    new_spans = []
    for (start, end) in spans:
        covered = set(xrange(start, end))
        if len(aligned_toks & covered) > 0:
            continue
        new_spans.append((start, end))
        aligned_toks |= covered
    return new_spans

# Here we want to pre-process AMR into a more generic form.
def categorized_amr(alignment_seq, amr, amr_statistics):

    aligned_set = set()
    tok_seq = amr.sent
    lemma_seq = amr.lems
    pos_seq = amr.poss

    concept_align_map = defaultdict(list)
    relation_align_map = defaultdict(list)

    # Initialize the alignment
    node_alignment, _ = alignment_utils.initializeAlignment(amr)

    node_to_toks, edge_to_toks, temp_aligned = alignment_utils.extractNodeMapping(
            alignment_seq, amr, concept_align_map, relation_align_map)

    temp_unaligned = set(xrange(len(pos_seq))) - temp_aligned

    aligned_toks = set()

    all_alignments = defaultdict(list)
    nodeid_to_frag = {}

    entity_toks = set()

    alignment_utils.alignEntities(tok_seq, amr, alignment_seq, nodeid_to_frag, entity_toks,
                                  aligned_toks, all_alignments, temp_unaligned, node_alignment)

    #Verbalization list
    verb_map = defaultdict(set)
    alignment_utils.alignVerbalization(tok_seq, lemma_seq, amr, VERB_LIST, all_alignments, verb_map,
                                       aligned_toks, node_alignment)

    aligned_nodes = set([node_idx for (node_idx, aligned) in enumerate(node_alignment) if aligned])

    alignment_utils.alignOtherConcepts(tok_seq, lemma_seq, amr, aligned_toks, aligned_nodes, node_to_toks,
                                       edge_to_toks, all_alignments)

    ##Based on the alignment from node index to spans in the string
    unaligned_set = set(xrange(len(pos_seq))) - aligned_toks
    unaligned_idxs = sorted(list(unaligned_set))
    logger.writeln("Unaligned tokens: %s" % (" ".join([tok_seq[i] for i in unaligned_idxs])))

    unaligned_nodes = amr.unaligned_nodes(aligned_nodes)
    logger.writeln("Unaligned vertices: %s" % " ".join([node.node_str() for node in unaligned_nodes]))

    assert len(tok_seq) == len(pos_seq)

    new_amr, _, span_to_type = AMRGraph.collapsedAMR(amr, all_alignments)

    collapsed_toks, collapsed_lem, collapsed_pos, new_alignment = collapseTokens(tok_seq, lemma_seq, pos_seq,
                                                                                     span_to_type, True)

    new_amr.set_sentence(collapsed_toks)
    new_amr.set_lemmas(collapsed_lem)
    new_amr.set_poss(collapsed_pos)

    return (new_amr, new_alignment)

