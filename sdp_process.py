#!/usr/bin/python
import sys
from collections import defaultdict
def processData(file, output_file=None):

    file_f = open(file)
    all_toks = []
    all_edge_maps = []
    all_piseqs = []

    toks = []
    lemmas = []
    poss = []
    edge_map = None
    pi_seq = []

    id_to_pred = {}
    preds = []

    num_preds = 0
    num_extra = 0

    for line in file_f:
        line = line.strip()

        if not line: #End of sentence, start a new sentence

            curr_edge_map = defaultdict(set)

            for pred_index in edge_map:
                curr_word_index = id_to_pred[pred_index]  #Get the word for the current predicate id
                for (rel, word_index) in edge_map[pred_index]:
                    curr_edge_map[curr_word_index].add(word_index)
                    curr_edge_map[word_index].add(curr_word_index)

            all_toks.append(toks)
            all_piseqs.append(pi_seq)
            for index in pi_seq:
                _ = curr_edge_map[index]
            all_edge_maps.append(curr_edge_map)

            toks = []
            lemmas = []
            pi_seq = []
            num_extra = 0
            edge_map = None
            num_preds = 0

            continue

        elif line[0] == '#': #Commented
            continue

        #Start processing each line
        splits = line.split()
        vertex_id = int(splits[0])
        pi_seq.append(vertex_id)
        toks.append(splits[1].strip())
        lemmas.append(splits[2].strip())
        poss.append(splits[3].strip())

        if splits[4].strip() == '+':
            is_top = True
        else:
            is_top = False

        if splits[5].strip() == '+':
            is_pred = True
            id_to_pred[num_preds] = int(splits[0])
            num_preds += 1
        else:
            is_pred = False

        #if splits[6].strip() != '_':
        #    print splits[6]

        temp_extra = len(splits) - 7
        if num_extra == 0:
            num_extra = temp_extra
        else:
            assert num_extra == temp_extra

        if edge_map is None:
            edge_map = defaultdict(list)

        for i in xrange(7, len(splits)):
            if splits[i].strip() != '_':
                edge_map[i-7].append((splits[i].strip(), int(splits[0])))

    if output_file != None:
        output_wf = open(output_file, 'w')
        n_sents = len(all_edge_maps)
        for i in xrange(n_sents):
            curr_edge_map = all_edge_maps[i]
            reprs = []
            for word_index in curr_edge_map:
                for tail_word_index in curr_edge_map[word_index]:
                    reprs.append("%d:%d" % (word_index, tail_word_index))
            print >> output_wf, (" ".join(reprs))
        output_wf.close()
    return all_toks, all_piseqs, all_edge_maps
    #for i in xrange(20):
    #    print all_edge_maps[i]
    #print num_extra
    #print num_preds

if __name__ == '__main__':
    processData(sys.argv[1], sys.argv[2])
