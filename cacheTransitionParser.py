#!/usr/bin/env python2.7

from collections import deque, defaultdict

from copy import copy, deepcopy

class State(object):
    def __init__(self, cache_, stack_, buffer_, connection_list_):
        self.cache = cache_
        self.stack = stack_
        self.buffer = buffer_
        self.connections = connection_list_

    def pop(self):
        assert len(self.stack) >= 2
        vertex = self.stack.pop()
        position = self.stack.pop()
        self.cache = self.cache[:position] + [vertex] + self.cache[position:-1]
        return position, vertex

    #From the cache, choose a vertex whose nearest connection in the buffer is farthest
    def push(self):

        def chooseVertex():
            push_position = 0
            push_vertex = 0
            max_dist = 0
            for (i, vertex) in enumerate(self.cache):
                if vertex not in self.connections: #Not connected to anything, infinity
                    push_position = i
                    push_vertex = vertex
                    break

                print self.cache, i, vertex
                #assert False, 'Should never get here'
                while self.connections[vertex] and self.connections[vertex][0] > len(self.buffer):
                    self.connections[vertex].popleft()

                if not self.connections[vertex]:
                    push_position = i
                    push_vertex = vertex
                    break

                if len(self.buffer) - self.connections[vertex][0] + 1 > max_dist:
                    max_dist = len(self.buffer) - self.connections[vertex][0] + 1
                    push_position = i
                    push_vertex = vertex


            return push_position, push_vertex

        push_position, push_vertex = chooseVertex()
        self.stack.append(push_position)
        self.stack.append(push_vertex)

        next_vertex = self.buffer.popleft()

        self.cache = self.cache[:push_position] + self.cache[push_position+1:] + [next_vertex]
        for (i, vertex) in enumerate(self.cache):
            if vertex == '$':
                continue
            while self.connections[vertex] and self.connections[vertex][0] > len(self.buffer):
                self.connections[vertex].popleft()
        return next_vertex

    @staticmethod
    def nextState(curr_state):
        next_state = deepcopy(curr_state)
        action = None
        if next_state.needsPop():
            next_state.pop()
            action = 'pop'
            introduced = None
        else:
            introduced = next_state.push()
            action = 'push'
        return action, next_state, introduced

    #Only pop when the rightmost vertex in the current cache
    #has no edge leading to vertices in the buffer
    def needsPop(self):
        r_vertex = self.cache[-1]
        if r_vertex == '$':
            return False

        #Pop when there is no vertex in the buffer that has edge connecting to the rightmost
        return len(self.connections[r_vertex]) == 0 or self.connections[r_vertex][0] > len(self.buffer)

    def finish(self):
        return len(self.buffer) == 0 and len(self.stack) == 0 and self.cache == (['$'] * len(self.cache))

    def printInfo(self):
        print 'Cache:', ' '.join([vertex if vertex == '$' else str(vertex) for vertex in self.cache])

    def info(self):
        return "%s %s" % ('Cache:', ' '.join([vertex if vertex == '$' else vertex for vertex in self.cache]))

class TreeNode(object):
    def __init__(self, state_):
        self.state = state_
        self.children = []
        self.parent = None
        self.introduced = None

    def add(self, node):
        self.children.append(node)

    def setParent(self, node):
        self.parent = node

    def setIntro(self, introduced_):
        self.introduced = introduced_

    def info(self):
        return self.state.info()

    def printInfo(self):
        self.state.printInfo()
        if self.introduced:
            print ("Introduce: %s" % self.introduced)


class TreeDecomposition(object):
    def __init__(self, root_):
        self.root = root_

    def width(self, node, edge_map):
        introduced_set = set()
        outside_set = set()
        stack = [node]

        counted_edges = set()
        while stack:
            curr_node = stack.pop()
            if curr_node.introduced is not None:
            #assert curr_node.introduced is not None
                introduced_set.add(curr_node.introduced)
            for child_node in curr_node.children:
                stack.append(child_node)

        for index in introduced_set:
            # needs to be verified if multiple subgraphs are disconnected.
            # try:
            #     if len(edge_map) > 1:
            #         assert index in edge_map
            # except:
            #     print index
            #     print edge_map
            #     print introduced_set
            #     sys.exit(1)
            for other_index in edge_map[index]:
                if other_index not in introduced_set:
                    outside_set.add(other_index)
                    counted_edges.add((index, other_index))
        #print outside_set
        return len(outside_set), counted_edges

    def computeTreeWidth(self, edge_map):
        max_width = 0
        stack = [self.root]
        while stack:
            curr_node = stack.pop()
            curr_width, _ = self.width(curr_node, edge_map)
            if curr_width > max_width:
                max_width = curr_width
            for child_node in curr_node.children:
                stack.append(child_node)
        return max_width

    def printTree(self, edge_map):
        stack = [(self.root, 0)]
        while stack:
            curr_node, depth = stack.pop()
            curr_str = ' ' * (depth*2)
            if curr_node.introduced is not None:
                curr_str += curr_node.introduced
            else:
                curr_str += "NONE"
            curr_width, counted_edges = self.width(curr_node, edge_map)
            node_info = curr_node.info()
            curr_str = "%s : %s : %d" % (curr_str, node_info, curr_width)
            if curr_width > 10:
                edges_str = ' '.join(['%s-%s' % (str(first), str(second)) for (first, second) in counted_edges])
                curr_str = "%s %s" % (curr_str, edges_str)
            print curr_str

            for child_node in reversed(curr_node.children):
                stack.append((child_node, depth+1))

    def printInfo(self):
        stack = [(self.root, 0)]
        while stack:
            curr_node, depth = stack.pop()
            print 'current depth: %d' % depth
            curr_node.printInfo()
            print 'Number of children: %d' % len(curr_node.children)
            for child_node in reversed(curr_node.children):
                stack.append((child_node, depth+1))

class DeterministicMachine(object):
    def __init__(self):
        return

    def start(self, pi_seq, edge_map_list, edge_map):

        introduced = set()

        n_vertices = len(pi_seq)
        init_cache = ['$'] * n_vertices
        pi_seq = deque(pi_seq)

        init_state = State(init_cache, [], pi_seq, edge_map_list)
        root_node = TreeNode(init_state)
        parent_node = root_node

        curr_state = init_state
        while not curr_state.finish():
            action, next_state, intro_vertex = State.nextState(curr_state)
            if action == 'push':
                new_node = TreeNode(next_state)
                new_node.setIntro(intro_vertex)
                parent_node.add(new_node)
                new_node.setParent(parent_node)
                parent_node = new_node
                introduced.add(intro_vertex)
            else: #Pop, return to the parent node
                parent_node = parent_node.parent
            curr_state = next_state

        td = TreeDecomposition(root_node)
        curr_width = td.computeTreeWidth(edge_map)
        # td.printInfo()

        return curr_width

