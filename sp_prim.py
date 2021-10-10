import sys
import math


# Prim's algorithm with Fibonacci heap

class Graph:
    def __init__(self, g):
        self.graph = g
        self.V = len(g)

    class Node:
        def __init__(self, index, key=sys.maxsize):
            self.source = None
            self.index = index
            self.key = key
            self.parent = self.child = self.left = self.right = None
            self.degree = 0
            self.mark = False

        def __repr__(self):
            if self.child is None:
                return f'Node({self.index})'
            else:
                return f'Node({self.index}): ' + self.child.__repr__()

        def add_child(self, node):
            self.degree += 1
            node.parent = self
            if self.child is None:
                self.child = node
                node.right = node.left = node
            else:
                self.child.left.right = node
                node.left = self.child.left
                node.right = self.child
                self.child.left = node

    class FibonacciHeap:
        def __init__(self):
            self.root = None
            self.minNode = None
            self.count = 0

        def __repr__(self):
            s = "Heap:\n"
            if self.root is not None:
                s += "->"
                s += str(self.root.key)
                p = self.root.right
                while p is not None and p != self.root:
                    s += "->"
                    s += str(p.key)
                    p = p.right
            return s

        def insert(self, node):
            self.count += 1
            self.insert_to_root(node)

        def extract_min(self):
            if self.count == 0:
                return None
            m = self.minNode
            # remove minNode from root list
            self.remove_from_root(m)
            self.minNode = None
            self.count -= 1
            # append child nodes to root list
            self.add_child_to_root(m)
            # consolidate the root list by combining nodes with same degree and update minNode
            self.consolidate()
            # return the old minNode
            return m

        def decrease_key(self, src, node, key):
            # decrease key
            if node.key <= key:
                return None
            node.key = key
            node.source = src
            p = node.parent
            if p is not None and p.key > node.key:
                self.cut(node)
                self.cascading_cut(p)
            if node.key < self.minNode.key:
                self.minNode = node

        def cut(self, node):
            # remove node from its parent's children
            if node == node.parent.child:
                if node.right == node:
                    node.parent.child = None
                else:
                    node.parent.child = node.right
            node.left.right = node.right
            node.right.left = node.left
            node.left = node.right = None
            # insert to root
            self.insert_to_root(node)
            # mark as false
            node.mark = False

        def cascading_cut(self, node):
            if node.parent is not None:
                if node.mark is False:
                    node.mark = True
                else:
                    self.cut(node)
                    self.cascading_cut(node.parent)

        def insert_to_root(self, node):
            if self.root is None:
                self.root = node
                node.right = node.left = node
                node.parent = None
                self.minNode = node
            else:
                if self.minNode.key > node.key:
                    self.minNode = node
                self.root.left.right = node
                node.left = self.root.left
                node.right = self.root
                self.root.left = node
                node.parent = None

        def remove_from_root(self, node):
            if node == self.root:
                if node.right == node:
                    self.root = None
                else:
                    self.root = node.right
            node.left.right = node.right
            node.right.left = node.left
            node.left = node.right = None

        def add_child_to_root(self, node):
            if node.child is not None:
                # set all children's parent to None
                temp = node.child
                while True:
                    temp.parent = None
                    temp = temp.right
                    if temp == node.child:
                        break
                # add to root
                if self.root is None:
                    self.root = node.child
                else:
                    self.root.left.right = node.child
                    node.child.left.right = self.root
                    node.child.left, self.root.left = self.root.left, node.child.left
                    node.degree = 0

        def consolidate(self):
            if self.root is None:
                return None
            arr = [None] * (math.ceil(2 * math.log2(self.count + 1)) + 1)
            while self.root is not None:
                temp = self.root
                self.remove_from_root(temp)
                while arr[temp.degree] is not None:
                    temp2 = arr[temp.degree]
                    arr[temp.degree] = None
                    if temp.key < temp2.key:
                        temp.add_child(temp2)
                    else:
                        temp2.add_child(temp)
                        temp = temp2
                arr[temp.degree] = temp
            for n in arr:
                if n is not None:
                    self.insert_to_root(n)

    def prim_mst_fheap(self):
        heap = self.FibonacciHeap()
        vertices = set()
        mst = set()

        for i in range(self.V):
            v = self.Node(i)
            vertices.add(v)
            heap.insert(v)

        heap.root.key = 0
        first_v = heap.extract_min()
        mst.add(first_v)
        vertices.remove(first_v)

        while heap.count != 0:
            for v1 in mst:
                for v2 in vertices:
                    heap.decrease_key(v1, v2, G[v1.index][v2.index])
            next_v = heap.extract_min()
            vertices.remove(next_v)
            mst.add(next_v)
            print(f'{next_v.source.index} - {next_v.index} (weight: {next_v.key})')

        return mst


G = [[0, 19, 5, 3, 2],
     [19, 0, 5, 9, 2],
     [5, 5, 0, 1, 6],
     [3, 9, 1, 0, 1],
     [2, 2, 6, 1, 0]]

g = Graph(G)
g.prim_mst_fheap()
