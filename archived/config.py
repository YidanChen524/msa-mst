"""
define global variables across files
"""


global score, gap, mapping
# define gap and score matrix
gap = 5
score = [[0, 5, 2, 5],
         [5, 0, 5, 2],
         [2, 5, 0, 5],
         [5, 2, 5, 0]]
# define how nucleotides are mapped to indices
mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3,
           'a': 0, 'c': 1, 'g': 2, 't': 3,
           'N': 0, 'R': 0, 'S': 1}
