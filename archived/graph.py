import sys


class Graph:

    def __init__(self, G):
        self.G = G
        self.V = len(self.G)

    def prim_list(self):

        mst = [False] * self.V
        src = [None] * self.V
        key = [sys.maxsize] * self.V

        # pick the first vertex
        mst[0] = True
        src[0] = -1
        key[0] = 0

        # find the minimum edge connecting to the current tree repeatedly
        src_v = 0
        for _ in range(self.V - 1):
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
            src_v = min_v

        # return the final tree, tree[i] = (src, dest, edge)
        tree = [None] * (self.V - 1)
        for i in range(1, self.V):
            tree[i - 1] = (src[i], i, key[i])
        return tree

    def prim_fibonacci_heap(self):

        class FibonacciHeap:

            class Node:

                def __init__(self, vertex, key=sys.maxsize, parent=None):
                    self.vertex = vertex
                    self.key = key
                    self.src = None
                    self.degree = 0
                    self.parent = None
                    self.children = None
                    self.right = self.left = None
                    self.mark = False

                def add_child(self, child):
                    # add a specific child node to self
                    self.degree += 1
                    child.parent = self
                    if self.children is None:
                        self.children = child
                        child.left = child
                        child.right = child
                    else:
                        self.children.left.right = child
                        child.left = self.children.left
                        child.right = self.children
                        self.children.left = child

                def remove_child(self, child):
                    # remove a specific child from self
                    self.degree -= 1
                    if self.degree == 0:
                        self.children = None
                    else:
                        self.children = child.right
                    child.parent = None
                    child.left.right = child.right
                    child.right.left = child.left
                    child.left = child.right = None
                    # return child
                    return child

            def __init__(self, G):
                self.G = G
                self.count = 0
                self.root = self.Node(-1)
                self.minNode = None

            def insert(self, vertex):
                # insert a new node into the vertices array
                node = self.Node(vertex)
                self.count += 1
                # add node to root
                self.root.add_child(node)
                # update minNode
                if self.minNode is None or node.key < self.minNode.key:
                    self.minNode = node

            def extract_min(self):
                # extract the minimum value from the heap
                m = self.root.remove_child(self.minNode)
                self.count -= 1
                self.minNode = None
                # add children of minNode to root
                while m.degree > 0:
                    child = m.remove_child(m.children)
                    self.root.add_child(child)
                # consolidate the root list by combining nodes with same degree and update minNode
                self.consolidate()
                # return (src, vertex, key) of the old minNode
                return m.src, m.vertex, m.key

            def consolidate(self):
                # merge nodes with same degree on the root list
                arr = [None] * self.root.degree
                while self.root.degree > 0:
                    n1 = self.root.remove_child(self.root.children)
                    n2 = arr[n1.degree]
                    while n2 is not None:
                        arr[n2.degree] = None
                        if n1.key < n2.key:
                            n1.add_child(n2)
                        else:
                            n2.add_child(n1)
                            n1 = n2
                        n2 = arr[n1.degree]
                    arr[n1.degree] = n1
                # add merged nodes back to root
                for n in arr:
                    if n is not None:
                        self.root.add_child(n)
                        # find the new minNode
                        if self.minNode is None or n.key < self.minNode.key:
                            self.minNode = n

            def decrease_key(self, src, node, key):
                # try decreasing the key of node to key
                # if success, set its src to src
                if node.key <= key or (key == 0 and src != -1):
                    return
                else:
                    node.key = key
                    node.src = src
                    parent = node.parent
                    if parent is not self.root and node.key < parent.key:
                        self.cut(node)
                        self.cascading_cut(parent)
                    if node.key < self.minNode.key:
                        self.minNode = node

            def cut(self, node):
                # cut the node from its parent, add to root, and mark it as False
                node.parent.remove_child(node)
                self.root.add_child(node)
                node.mark = False

            def cascading_cut(self, node):
                # if node's parent is marked as True, cut it and cascading cut its grandparent
                # else, mark it as True
                parent = node.parent
                if parent is not self.root:
                    if node.mark is False:
                        node.mark = True
                    else:
                        self.cut(node)
                        self.cascading_cut(parent)

            def initialize(self):
                # pick the first node and set its key to 0
                self.decrease_key(-1, self.root.children, 0)

            def decrease_all_connecting_keys(self, src):
                # decrease keys of all nodes connecting to src
                self.nodes = [None] * self.count
                self.index = 0
                self.traverse_heap(self.root)
                for node in self.nodes:
                    self.decrease_key(src, node, self.G[src][node.vertex])

            def traverse_heap(self, node):
                # traverse the heap and store current nodes into self.nodes
                if node is not None:
                    if node is not self.root:
                        self.nodes[self.index] = node
                        self.index += 1
                    # traverse horizontally
                    child = node.children
                    while child is not None:
                        # traverse vertically
                        self.traverse_heap(child)
                        child = child.right
                        if child is node.children:
                            break

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


G1 = [[0, 19, 5, 3, 2],
      [19, 0, 5, 9, 2],
      [5, 5, 0, 1, 6],
      [3, 9, 1, 0, 1],
      [2, 2, 6, 1, 0]]

G2 = [[0, 2, 0, 6, 0],
      [2, 0, 3, 8, 5],
      [0, 3, 0, 0, 7],
      [6, 8, 0, 0, 9],
      [0, 5, 7, 9, 0]]

if __name__ == "__main__":
    g = Graph(G2)
    print(g.prim_list())
    print(g.prim_fibonacci_heap())
