"""
Simulating fasta files with different number of sequences and length
"""

import random

DNA = 'ACGT'
LINE_WIDTH = 60


def simulate_string(m):
    """Simulate a DNA sequence of length m"""
    nucleotides = [random.choice(DNA) for _ in range(m)]
    lines = []
    for i in range(0, m, LINE_WIDTH):
        lines.append(''.join(nucleotides[i:(i + LINE_WIDTH)]))
    return '\n'.join(lines)


def simulate_sequences(n, m):
    """Simulate n sequences of length m."""
    with open(f"test_{n}_{m}.fa", "w") as f:
        for i in range(n):
            f.write(f">seq{i}\n")
            f.write(simulate_string(m))
            f.write("\n\n")


if __name__ == '__main__':
    for i in range(10, 101, 10):
        for j in range(10, 101, 10):
            simulate_sequences(i, j)
