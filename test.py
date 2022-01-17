from seqs.seqs import Seqs

s = Seqs("tests/random_seqs/test_10_10.fa")
s.global_align("approx")
s.output_alignments()
