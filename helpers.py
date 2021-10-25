"""
helper functions
"""


def parse_fasta(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    names = []
    seqs = []
    for i in range(len(lines)):
        if lines[i][0] == '>':
            names.append(lines[i][1:].strip('\n'))
            seq = ""
            while i < len(lines) - 1 and lines[i + 1][0] != '>':
                seq += lines[i + 1].strip('\n')
                i += 1
            seqs.append(seq)
    return names, seqs


def string_concat(*s):
    """concatenate strings in O(n)"""
    return "".join(s)
