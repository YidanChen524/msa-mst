import sys
from .fheap import FibonacciHeap


class Graph:

    def __init__(self, G):
        self.G = G
        self.V = len(self.G)

    def prim_list_mst(self):
        mst = [False] * self.V
        src = [None] * self.V
        key = [sys.maxsize] * self.V
        # pick the first vertex
        mst[0] = True
        src[0] = -1
        key[0] = 0
        # find the minimum edge connecting to the current tree repeatedly
        src_v = 0
        tree = [None] * (self.V - 1)
        for i in range(self.V - 1):
            # update keys related to the previously picked vertices
            for v in range(self.V):
                if mst[v] is False and 0 < self.G[src_v][v] < key[v]:
                    key[v] = self.G[src_v][v]
                    src[v] = src_v
            # find the new vertex with minimum key
            min_key = sys.maxsize
            min_v = None
            for v in range(self.V):
                if key[v] < min_key and mst[v] is False:
                    min_key = key[v]
                    min_v = v
            mst[min_v] = True
            tree[i] = (src_v, min_v, min_key)
            src_v = min_v
        # return the final tree, tree[i] = (src, dest, edge)
        return tree

    def prim_fheap_mst(self):
        # initialize heap and insert all vertices into the heap
        heap = FibonacciHeap(self.G)
        for i in range(self.V):
            heap.insert(i)
        # pick the first vertex
        heap.initialize()
        heap.extract_min()
        heap.decrease_all_connecting_keys(0)
        # extend the minimum tree repeatedly
        tree = [None] * (self.V - 1)
        for i in range(1, self.V):
            src, vertex, key = heap.extract_min()
            heap.decrease_all_connecting_keys(vertex)
            tree[i-1] = (src, vertex, key)
        # return the mst tree
        return tree
