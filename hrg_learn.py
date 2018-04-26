#!/usr/bin/python
import sys
import os
from amr_utils import *
from HRGSample import *
import argparse
from date_extraction import *
from treewidth import depthFirstPi, buildPiSeq, buildRightNeighbors, ordered_edges
from amrtreewidth import instanceTreewidth
from categorize_amr import categorized_amr
from hgraph import *
from ioutil import loadDependency
import alignment_utils

def mergeSpans(index_to_spans):
    new_idx_to_spans = {}
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
        new_idx_to_spans[index] = new_span_list
    return new_idx_to_spans

def getDateAttr(frag):
    date_relations = set(['time', 'year', 'month', 'day', 'weekday', 'century', 'era', 'decade', 'dayperiod', 'season', 'timezone'])
    assert len(frag.roots) == 1
    root_idx = frag.roots[0]
    amr_graph = frag.graph
    root_node = amr_graph.nodes[root_idx]

    index_to_attr = {}

    attr_indices = set()
    for edge_idx in root_node.v_edges:
        curr_edge = amr_graph.edges[edge_idx]
        if curr_edge.label in date_relations:
            attr_indices.add(curr_edge.tail)
            index_to_attr[curr_edge.tail] = curr_edge.label.upper()
    return (attr_indices, index_to_attr)


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

        dump_file(pred_f, self.num_predicates)
        dump_file(non_pred_f, self.num_nonpredicate_vals)
        dump_file(const_f, self.num_consts)
        dump_file(entity_f, self.num_entities)
        dump_file(named_entity_f, self.num_named_entities)

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

def is_projective(n, edge_map):
    right_map, right_most_map = buildRightNeighbors(n, edge_map)
    for i in xrange(n):
        if i in right_map:
            for j in right_map[i]:
                assert j > i
                for idx in xrange(i+1, j):
                    if idx in right_map:
                        assert idx in right_most_map
                        if right_most_map[idx] > j:
                            return False
    return True

# extract all the binarized fragments combinations for AMR
# each chart item is a set of fragments are consistent with a span of aligned strings
# Each chart item is a unique subgraph.
def optimal_tree_decomposition(n, edge_map, edge_label_map,
                              inside_rule_dist, inorder_rule_dist,
                              nolabel=False, nodir=False):

    def print_decomposition(bp, wd, bgs):
        stack = [(0, n, 0)]
        while stack:
            l, r, level = stack.pop()
            curr_repr = "%s(%d %d width: %d %s)" % ("  "*level, l, r, wd[l][r], str(bgs[(l, r)]))
            print curr_repr
            assert (l, r) in bp
            if bp[(l, r)] is not None:
                for (new_l, new_r) in reversed(bp[(l, r)]):
                    stack.append((new_l, new_r, level+1))

    # Given a span, compute the vertices that connect with the outside
    def compute_ext_set(span_set, edge_map):
        dockings = set()
        for vertex in span_set:
            if vertex in edge_map:
                outside_set = edge_map[vertex] - span_set
                if len(outside_set) > 0:
                    dockings.add(vertex)
        return dockings

    # Given a span, and the edge maps, compute the outside vertex set efficiently.
    def compute_outside_set(span_set, dockings, edge_map):
        connection_set = set()
        for vertex in dockings:
            assert vertex in edge_map
            connection_set |= edge_map[vertex]
        return connection_set - span_set

    def inside_outside(start, end, edge_map):
        span_set = set(xrange(start, end))
        dockings = compute_ext_set(span_set, edge_map)
        outside_set = compute_outside_set(span_set, dockings, edge_map)
        return dockings, outside_set

    def extract_derivation(bp, bgs, label_map):
        """
        param bp: backpointers for tracking the children information.
        param bgs: the vertices contained in each bag, ordered list.
        param labeled: whether each edge is labeled or not.
        return: the list of derivation rules.
        """
        stack = [(0, n, 0, [])]
        derivation = []
        visited = set()
        while stack:
            l, r, level, ext_list = stack.pop()
            curr_bag = bgs[(l, r)]
            curr_hg = Hypergraph()
            curr_hg.init_vertices(curr_bag)
            vertex_lab_map = curr_hg.vertex_idx
            ext_labs = []

            # Map vertex order index to current label space.
            for vertex in ext_list:
                ext_labs.append(vertex_lab_map[vertex])

            curr_hg.setExtList(ext_labs)

            if bp[(l, r)] is not None:
                next_lists = []
                # Deal with all the nonterminal edges.
                for (child_l, child_r) in bp[(l, r)]:
                    child_bag = bgs[(child_l, child_r)]
                    child_connect = curr_bag & child_bag
                    try:
                        assert len(child_connect) > 0
                    except:
                        print "Weird decomposition here"
                        # print_decomposition(bp, width, bgs)
                        # return derivation

                    connect_list = sorted([vertex for vertex in child_connect])
                    next_lists.append((child_l, child_r, level+1, connect_list))
                    child_lab = "X%d" % len(connect_list)
                    connect_labs = [vertex_lab_map[vertex] for vertex in connect_list]
                    child_edge = Hyperedge(connect_labs, child_lab)
                    curr_hg.addEdge(child_edge)
                stack.extend(next_lists[::-1])
            else:
                assert r == l + 1 or r == -1
                # curr_lab = "X-%s" % vertex_labs[l]
                curr_lab = "X1"
                curr_edge = Hyperedge([vertex_lab_map[l]], curr_lab)
                curr_hg.addEdge(curr_edge)

                if len(curr_bag) == 1:
                    continue

            # Then we process the terminal edges.
            # These edges start and end with vertices in the current bag.
            for start in curr_bag:
                for end in curr_bag:
                    if (start, end) in label_map and ((start, end) not in visited):
                        if start >= end:
                            continue
                        # assert start < end
                        if nolabel:
                            if nodir:
                                curr_lab = "E0"
                            else:
                                curr_lab = label_map[(start, end)][0]
                        else:
                            curr_lab = label_map[(start, end)]
                        start_lab = vertex_lab_map[start]
                        end_lab = vertex_lab_map[end]
                        curr_edge = Hyperedge([start_lab, end_lab], curr_lab)
                        curr_hg.addEdge(curr_edge)
                        visited.add((start, end))

            curr_repr = (level, str(curr_hg))
            derivation.append(curr_repr)
        return derivation

    # projective = is_projective(n, edge_map)
    # if not projective:
    #     print "Non-projective dependency tree"
    #     print edge_map

    chart = [[None for _ in range(n+1)] for _ in range(n+1)]
    inside_width = [[n for _ in range(n+1)] for _ in range(n+1)]
    inside_bags = {}
    inside_pointers = {}

    inorder_width = [[n for _ in range(n+1)] for _ in range(n+1)]
    for i in range(n+1):
        inorder_width[i][i] = 0
    inorder_bags = {}
    inorder_pointers = {}

    # The leaves of the forest are ordered vertices of the graph.
    for i in xrange(n):
        j = i + 1

        chart[i][j] = inside_outside(i, j, edge_map)

        # Inside tree decomposition.
        inside_width[i][j] = 1
        inside_bags[(i, j)] = set([i])
        inside_pointers[(i, j)] = None

        outside_set = chart[i][j][1]
        out_bag = outside_set | set([i])

        # Initialize inorder tree decomposition.
        # Anchored leaf node.
        inorder_width[i][j] = len(outside_set) + 1
        inorder_bags[(i, j)] = out_bag
        inorder_pointers[(i, j)] = None

    for span in xrange(2, n+1):
        for i in xrange(0, n):
            j = i + span
            if j > n:
                continue

            chart[i][j] = inside_outside(i, j, edge_map)
            curr_ext_set, curr_outside_set = chart[i][j]

            # First compute the treewidth of inside tree decomposition.
            for k in xrange(i+1, j):
                assert chart[i][k] is not None and chart[k][j] is not None

                # Then we update the width stats
                bsize = len(chart[i][k][0] | chart[k][j][0])
                # assert bsize > 0
                curr_width = max(inside_width[i][k], inside_width[k][j], bsize)
                if curr_width <= inside_width[i][j]:
                    inside_width[i][j] = curr_width
                    inside_pointers[(i, j)] = [(i, k), (k, j)]
                    inside_bags[(i, j)] = chart[i][k][0] | chart[k][j][0]

            # Here we compute the treewidth of in-order tree decomposition.
            # First unanchored binary nodes.
            for k in range(i, j):
                left_set = set(range(i, k))
                right_set = set(range(k, j))
                if k == i:
                    left_sect = set()
                else:
                    left_sect = chart[i][k][1] & right_set
                right_sect = chart[k][j][1] & left_set

                if k > i and len(left_sect) == 0 and len(right_sect) == 0:
                    curr_width = max(len(curr_outside_set), inorder_width[i][k],
                                     inorder_width[k][j])
                    if curr_width <= inorder_width[i][j]:
                        inorder_width[i][j] = curr_width
                        inorder_pointers[(i, j)] = [(i, k), (k, j)]
                        inorder_bags[(i, j)] = curr_outside_set

                right_set = set(range(k+1, j))
                if i < k:
                    left_sect = chart[i][k][1] & right_set
                else:
                    left_sect = set()

                if k < j-1:
                    right_sect = chart[k+1][j][1] & left_set
                else:
                    right_sect = set()
                if len(left_sect) == 0 and len(right_sect) == 0:
                    curr_width = max(len(curr_outside_set)+1, inorder_width[i][k],
                                     inorder_width[k+1][j])

                    if curr_width <= inorder_width[i][j]:
                        inorder_width[i][j] = curr_width
                        inorder_bags[(i, j)] = curr_outside_set | set([k])
                        if k == i:
                            inorder_pointers[(i, j)] = [(i, -1), (i+1, j)]
                            inorder_pointers[(i, -1)] = None
                            inorder_bags[(i, -1)] = set([i])
                        elif k == j-1:
                            inorder_pointers[(i, j)] = [(i, j-1), (j-1, -1)]
                            inorder_pointers[(j-1, -1)] = None
                            inorder_bags[(j-1, -1)] = set([j-1])
                        else:
                            inorder_pointers[(i, j)] = [(i, k), (k, -1), (k+1, j)]
                            inorder_pointers[(k, -1)] = None
                            inorder_bags[(k, -1)] = set([k])

    assert chart[0][n] is not None, "None found for current sentence."

    # print "Inside bag information"
    # print inside_bags
    curr_der = extract_derivation(inside_pointers, inside_bags, edge_label_map)
    curr_der_repr = "%s\n" % "\n".join([("  "*level+ rule_str) for
                                        (level, rule_str) in curr_der])
    for (_, rule_str) in curr_der:
        inside_rule_dist[rule_str] += 1

    if n < 10:
        print "Inside tree decomposition details:"
        print_decomposition(inside_pointers, inside_width, inside_bags)

    print >> inside_rule_f, curr_der_repr

    curr_der = extract_derivation(inorder_pointers, inorder_bags, edge_label_map)
    curr_der_repr = "%s\n" % "\n".join([("  "*level+ rule_str) for
                                        (level, rule_str) in curr_der])
    for (_, rule_str) in curr_der:
        inorder_rule_dist[rule_str] += 1

    print >> inorder_rule_f, curr_der_repr

    if n < 10:
        print "Inorder tree decomposition details:"
        print_decomposition(inorder_pointers, inorder_width, inorder_bags)

    return inside_width[0][n]-1, inorder_width[0][n]-1


def construct_forest(args):
    sys.setrecursionlimit(sys.getrecursionlimit() * 30)

    # Rule stats files for inside and inorder tree decomposition.
    global inside_rule_f
    global inorder_rule_f

    os.system("mkdir -p %s" % args.output_dir)

    if args.nolabel:
        if args.nodir:
            inside_rule_f = open(os.path.join(args.output_dir, "inside_derivation_nolabel_nodir.txt"), "w")
            inorder_rule_f = open(os.path.join(args.output_dir, "inorder_derivation_nolabel_nodir.txt"), "w")
            inside_stat_f = open(os.path.join(args.output_dir, "inside_rule_stats_nolabel_nodir.txt"), 'w')
            inorder_stat_f = open(os.path.join(args.output_dir, "inorder_rule_stats_nolabel_nodir.txt"), 'w')
        else:
            inside_rule_f = open(os.path.join(args.output_dir, "inside_derivation_dironly.txt"), "w")
            inorder_rule_f = open(os.path.join(args.output_dir, "inorder_derivation_dironly.txt"), "w")
            inside_stat_f = open(os.path.join(args.output_dir, "inside_rule_stats_dironly.txt"), 'w')
            inorder_stat_f = open(os.path.join(args.output_dir, "inorder_rule_stats_dironly.txt"), 'w')
    else:
        inside_rule_f = open(os.path.join(args.output_dir, "inside_derivation_labeled.txt"), "w")
        inorder_rule_f = open(os.path.join(args.output_dir, "inorder_derivation_labeled.txt"), "w")
        inside_stat_f = open(os.path.join(args.output_dir, "inside_rule_stats_labeled.txt"), 'w')
        inorder_stat_f = open(os.path.join(args.output_dir, "inorder_rule_stats_labeled.txt"), 'w')

    # Global width statistics.
    n_success = 0.0
    total_width = 0.0
    max_width = 0.0
    width_dist = defaultdict(int)
    inside_rule_dist = defaultdict(int)

    total_inorder_width = 0.0
    max_inorder_width = 0.0
    inorder_dist = defaultdict(int)
    inorder_rule_dist = defaultdict(int)

    total_cache_width = 0.0
    max_cache_width = 0.0
    cache_dist = defaultdict(int)

    num_projective = 0.0

    amr_format = (args.format.lower() == "amr")

    if amr_format:
        amr_file = os.path.join(args.data_dir, 'amr')
        alignment_file = os.path.join(args.data_dir, 'alignment')
        lemma_file = os.path.join(args.data_dir, 'lemmatized_token')
        tok_file = os.path.join(args.data_dir, 'token')
        pos_file = os.path.join(args.data_dir, 'pos')

        amr_graphs = load_amr_graphs(amr_file)
        alignments = [line.strip().split() for line in open(alignment_file, 'r')]
        toks = [line.strip().split() for line in open(tok_file, 'r')]
        lemmas = [line.strip().split() for line in open(lemma_file, 'r')]
        poss = [line.strip().split() for line in open(pos_file, 'r')]

        assert len(amr_graphs) == len(alignments) and len(amr_graphs) == len(toks) and len(amr_graphs) == len(poss), (
            '%d %d %d %d %d' % (len(amr_graphs), len(alignments), len(toks), len(poss)))

        num_self_cycle = 0

        amr_statistics = AMR_stats()

        if args.use_stats:
            amr_statistics.loadFromDir(args.stats_dir)
            print amr_statistics
        else:
            os.system('mkdir -p %s' % args.stats_dir)
            amr_statistics.collect_stats(amr_graphs)
            amr_statistics.dump2dir(args.stats_dir)
    else: # In conll dependency format.
        toks, all_deps = loadDependency(args.conll_file)

    random.seed(0)
    for (sent_idx, tok_seq) in enumerate(toks):

        print >> inside_rule_f, ("Sentence %d" % sent_idx)
        print >> inside_rule_f, ("Original sentence: %s" % " ".join(tok_seq))

        print >> inorder_rule_f, ("Sentence %d" % sent_idx)
        print >> inorder_rule_f, ("Original sentence: %s" % " ".join(tok_seq))

        # if sent_idx != 4187:
        #     continue

        print 'Sentence #%d' % sent_idx
        print "Original sentence:", " ".join(tok_seq)

        if amr_format:
            amr = amr_graphs[sent_idx]
            print >> inside_rule_f, ("AMR: \n%s" % str(amr))
            print >> inorder_rule_f, ("AMR: \n%s" % str(amr))
            lemma_seq, pos_seq, alignment_seq = lemmas[sent_idx], poss[sent_idx], alignments[sent_idx]

            amr.setStats(amr_statistics)

            if amr.check_self_cycle():
                num_self_cycle += 1

            amr.set_sentence(tok_seq)
            amr.set_lemmas(lemma_seq)
            amr.set_poss(pos_seq)

            new_amr, new_alignment = categorized_amr(alignment_seq, amr, amr_statistics)

            span_idxes = []
            for node_idx in new_alignment:
                for (start, end, _) in new_alignment[node_idx]:
                    span_idxes.append((start, end, node_idx))
            sorted_idxes = sorted(span_idxes, key=lambda x: (x[0], x[1]))

            sorted_idxes = [z for (x, y, z) in sorted_idxes]

            # The pi seq is ordered, while edges are already transformed.
            pi_seq, edge_map, edge_label_map = buildPiSeq(new_amr, sorted_idxes)

            print "categorized sentence:", " ".join(new_amr.sent)
            print str(new_amr)

            print pi_seq
            num_vertices = len(pi_seq)

            # Write down the sequence of vertices.
            vertex_labels = [new_amr.nodes[node_id].node_label() for node_id in pi_seq]
            print "Vertex order:", " ".join(["%d:%s" % (id, l) for (id, l) in enumerate(vertex_labels)])
            print edge_label_map

        else:
            num_vertices, edge_map, edge_label_map = all_deps[sent_idx]
            print num_vertices, edge_map, edge_label_map

        projective = is_projective(num_vertices, edge_map)
        if not projective:
            print "Non-projective dependency tree"
            print edge_map
        else:
            num_projective += 1.0

        curr_width, inorder_width = optimal_tree_decomposition(num_vertices, edge_map, edge_label_map, inside_rule_dist,
                                                               inorder_rule_dist, args.nolabel, args.nodir)
        pi_seq = list(range(num_vertices))
        cache_width = instanceTreewidth(pi_seq, edge_map)

        if curr_width is not None:
            if len(pi_seq) < 7:

                if inorder_width < curr_width:
                    print "inorder better than inside"
                elif inorder_width > curr_width:
                    print "inside better than inorder"
            print 'current cache width:', cache_width
            print 'current inside width:', curr_width
            print 'current inorder width:', inorder_width

            if projective:
                assert curr_width <= 1, ("A projective tree with inside width of %d" % curr_width)
                assert inorder_width <=1, ("A projective tree with in-order width of %d" % inorder_width)

            print "\n"

            n_success += 1.0

            total_width += curr_width
            if curr_width > max_width:
                max_width = curr_width
            width_dist[curr_width] += 1

            total_inorder_width += inorder_width
            if inorder_width > max_inorder_width:
                max_inorder_width = inorder_width
            inorder_dist[inorder_width] += 1

            total_cache_width += cache_width
            if cache_width > max_cache_width:
                max_cache_width = cache_width
            cache_dist[cache_width] += 1

        else:
            print 'Failed to construct forest'

    sorted_inside_rules = sorted(inside_rule_dist.items(), key=lambda x: -x[1])
    for (rule, count) in sorted_inside_rules:
        print >> inside_stat_f, "%s\t\t%d" % (rule, count)
    inside_stat_f.close()

    sorted_inorder_rules = sorted(inorder_rule_dist.items(), key=lambda x: -x[1])
    for (rule, count) in sorted_inorder_rules:
        print >> inorder_stat_f, "%s\t\t%d" % (rule, count)
    inorder_stat_f.close()

    print 'Analyzed a total of %d sentences.' % int(n_success)
    print 'Projective ratio is:', num_projective/n_success
    print "cache transition treewidth:"
    sorted_cache_width = sorted(cache_dist.items(), key=lambda x: -x[1])
    for (width, count) in sorted_cache_width:
        print 'width: %d, count: %d' % (width, count)
    print 'Average cache width:', total_cache_width/n_success

    print "inside treewidth:"
    sorted_width = sorted(width_dist.items(), key=lambda x: -x[1])
    for (width, count) in sorted_width:
        print 'width: %d, count: %d' % (width, count)
    print 'Average inside width:\n', total_width/n_success

    print "inorder treewidth:"
    sorted_inorder_width = sorted(inorder_dist.items(), key=lambda x: -x[1])
    for (width, count) in sorted_inorder_width:
        print 'width: %d, count: %d' % (width, count)
    print 'Average inorder treewidth:\n', total_inorder_width/n_success

    inside_rule_f.close()
    inorder_rule_f.close()

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()

    argparser.add_argument("--use_stats", action="store_true", help="if to use stats")
    argparser.add_argument("--refine", action="store_true", help="if to refine the nonterminals")
    argparser.add_argument("--nolabel", action="store_true", help="whether to ignore edge label")
    argparser.add_argument("--nodir", action="store_true", help="whether to ignore edge direction")
    argparser.add_argument("--random", action="store_true", help="if use random order")
    argparser.add_argument("--reversed", action="store_true", help="if to reverse the pi sequence")
    argparser.add_argument("--depth", action="store_true", help="if to reverse the pi sequence")
    argparser.add_argument("--data_dir", type=str, help="the data directory for dumped AMR graph "
                                                        "objects, alignment and tokenized sentences")
    argparser.add_argument("--conll_file", type=str, help="the conll format dependency file")
    argparser.add_argument("--output_dir", type=str, help="the output directory for saving the constructed forest")
    argparser.add_argument("--format", type=str, help="the format for the input file")
    argparser.add_argument("--stats_dir", type=str, help="the statistics directory")
    argparser.add_argument("--nsplit", type=int, help="if to split the forest into multiple")

    args = argparser.parse_args()
    construct_forest(args)
