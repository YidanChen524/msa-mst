"""
Simulating fasta files with different number of sequences and length for tests
"""

import random
import os

DNA = 'ACGT'
LINE_WIDTH = 60


def simulate_random_string(m):
    """Simulate a DNA sequence of length m randomly"""
    nucleotides = [random.choice(DNA) for _ in range(m)]
    lines = []
    for i in range(0, m, LINE_WIDTH):
        lines.append(''.join(nucleotides[i:(i + LINE_WIDTH)]))
    return '\n'.join(lines)


def simulate_related_string(template, proportion):
    """Simulate a DNA sequence given a template sequence and maximum edit proportion"""
    m = len(template)
    edits = int(m*proportion)
    nucleotides = list(template)
    for _ in range(edits):
        pos = random.randrange(m)
        nucleotides[pos] = random.choice(DNA)
    lines = []
    for i in range(0, m, LINE_WIDTH):
        lines.append(''.join(nucleotides[i:(i + LINE_WIDTH)]))
    return '\n'.join(lines)


def simulate_random_sequences(n, m):
    """Simulate n sequences of length m."""
    with open(f"{os.path.dirname(os.path.abspath(__file__))}/random_seqs/test_{n}_{m}.fa", "w") as f:
        for i in range(n):
            f.write(f">seq{i}\n")
            f.write(simulate_random_string(m))
            f.write("\n\n")


def simulate_related_sequences(n, m):
    """Simulate n sequences of length m."""
    with open(f"{os.path.dirname(os.path.abspath(__file__))}/related_seqs/test_{n}_{m}.fa", "w") as f:
        template_seq = simulate_random_string(m)
        f.write(f">seq0\n{template_seq}\n\n")
        for i in range(1, n):
            f.write(f">seq{i}\n")
            f.write(simulate_related_string(template_seq, i/n))
            f.write("\n\n")


if __name__ == '__main__':
    for i in range(10, 301, 10):
        for j in range(10, 101, 10):
            # simulate_random_sequences(i, j)
            simulate_related_sequences(i, j)
