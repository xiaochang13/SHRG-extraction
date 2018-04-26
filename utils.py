#!/usr/bin/python
import re
import sys
from collections import defaultdict
symbols = set("'\".-!#&*|\\/")
re_symbols = re.compile("['\".\'`\-!#&*|\\/@=\[\]]")
def allSymbols(s):
    return re.sub(re_symbols, "", s) == ""

def equals(s1, s2):
    s2_sub = s2.replace("-RRB-", ")").replace("-LRB-", "(")
    s2_sub = re.sub(re_symbols, "", s2_sub)
    s1 = s1.replace("\xc2\xa3", "#")
    s1_sub = re.sub(re_symbols, "", s1)

    return s1_sub.lower() == s2_sub.lower()
    #     return True

    # if len(s1) > 6 and len(s2) > 6:
    #     return s1[:6] == s2[:6] or s1[-6:] == s2[-6:]

    # return False

#Search s in seq from index
def searchSeq(s, seq, index):
    for i in xrange(index, len(seq)):
        if s == "@-@":
            if "".join(seq[i:i+3]) == "@-@":
                return (i, i+3)
        if s == "\"":
            if "".join(seq[i:i+2]) == "\'\'":
                return (i, i+2)
            if seq[i] == "\'\'":
                return (i, i+1)
        if allSymbols(s):
            if seq[i] == s:
                return (i, i+1)
        if allSymbols(seq[i]):
            continue
        for j in xrange(i+1, len(seq)+1):
            if allSymbols(seq[j-1]):
                continue
            curr_seq = "".join(seq[i:j])
            if equals(s, curr_seq):
                return (i, j)
            #else:
            #    print "test:", s, curr_seq, re.sub(re_symbols, "", curr_seq), "False"
    #print s, re.sub(re_symbols, "", s)
    #print " ".join(seq[index:])
    return (-1, -1)

def loadTokens(input_file):
    ret = []
    with open(input_file, 'r') as input_f:
        for line in input_f:
            ret.append(line.strip().split())
        input_f.close()
    return ret

#Based on some heuristic tokenize rules, map the original toks to the tokenized result
#Also map the toks to the new alignment toks
def mergeToks(orig_toks, tokenized_seq, all_alignments, sent_index):
    new_alignment = defaultdict(list)
    matched_index = 0
    triple_list = []
    for index in all_alignments:
        for (start, end, wiki_label) in all_alignments[index]:
            triple_list.append((index, (start, end, wiki_label)))
    sorted_alignments = sorted(triple_list, key=lambda x: (x[1][0], x[1][1]))

    visited = set()
    for (index, (start, end, wiki_label)) in sorted_alignments:
        for i in xrange(start, end):
            if i in visited:
                continue
            if i < end:
                curr_span = "".join(orig_toks[i: end])
                if allSymbols(curr_span):
                    print curr_span, wiki_label
                    break

                (new_start, new_end) = searchSeq(curr_span, tokenized_seq, matched_index)
                if new_start == -1:
                    print ("Something wrong here in %d" % sent_index)
                    print curr_span
                    print orig_toks
                    print tokenized_seq
                    print matched_index, tokenized_seq[matched_index]

                    sys.exit(1)
                visited |= set(xrange(i, end))
                matched_index = new_end
                new_alignment[index].append((new_start, new_end, wiki_label))
                #print " ".join(orig_toks[i:end]), "  :  ", " ".join(tokenized_seq[new_start:new_end])
    return new_alignment

