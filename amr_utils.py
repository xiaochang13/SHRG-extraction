#!/usr/bin/python
import cPickle
import sys
import os
import amr_graph
from amr_graph import *
from re_utils import *

def get_amr_line(input_f):
    """Read the amr file. AMRs are separated by a blank line."""
    cur_amr=[]
    has_content=False
    for line in input_f:
      if line[0]=="(" and len(cur_amr)!=0:
         cur_amr=[]
      if line.strip()=="":
         if not has_content:
            continue
         else:
            break
      elif line.strip().startswith("#"):
        # omit the comment in the AMR file
        continue
      else:
         has_content=True
         cur_amr.append(delete_pattern(line.strip(), '~e\.[0-9]+(,[0-9]+)*'))
         #cur_amr.append(line.strip())
    return "".join(cur_amr)


# Load a list of amr graph objects
def load_amr_graphs(amr_file):
    f = open(amr_file, 'r')
    amr_line = get_amr_line(f)
    graphs = []
    while amr_line and amr_line != '':
        amr_graph = AMRGraph(amr_line)
        graphs.append(amr_graph)
        amr_line = get_amr_line(f)

    return graphs

if __name__ == '__main__':
    s = '(f / foolish :condition (d / do-02  :ARG0 i) :domain (i / i))'

