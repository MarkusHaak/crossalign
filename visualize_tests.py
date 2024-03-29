import unittest
import os
import pickle
import pandas as pd
from crossalign import *

def write_fastq(fn, seq_id, seq, mode="w"):
    with open(fn, mode) as f:
        qual_str = "!"*len(seq)
        print(f"@{seq_id}\n{seq}\n+\n{qual_str}", file=f)

def write_fasta(fn, seq_id, seq, mode="w"):
    with open(fn, mode) as f:
        print(f">{seq_id}\n{seq}", end="", file=f)

class SimpleAlignmentVisualizationCase(unittest.TestCase):

    def setUp(self):
        # for these test, no external files are needed
        global args
        args = parse_args(
            [
                '--reads', 'None',
                '--genome', 'None',
                '--adapter', 'None',
                '--mean', '1.15', # normally autom. determined
                '--std', '0.1',   # normally autom. determined
                #'--mismatch', '-2',
                #'--gap_open', '0',
                #'--gap_extension', '-1',
            ])
        init_ctypes_arrays(args)
        ct.CDLL(clib).init(args.sequence_type == 'nucl')

    def test_align_1(self):
        """high ambiguity"""
        free_gap = [True, True]
        query_seq = "GATGATCC"
        ref1      = "GAGACGG"
        ref2      = "CCCGACC"
    
        score = c_align(args, query_seq, ref1, ref2, free_gap)
        transitions_list = retrieve_transitions_list(query_seq, ref1, ref2, '+', '+', 0, 0, 0, 0)
    
        for i,(ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2) in enumerate(transitions_list):
            cigar = "".join([inflate_cigar(cg_ref1), inflate_cigar(cg_gap), inflate_cigar(cg_ref2)])
            print(f"\n#{i+1}: transition {query_ts}->{query_te}, score: {score}, cigar: {cigar}")
            print_crossalignment(query_seq, ref1, ref2, cigar, fna_ref1, fa_ref2, query_ts, query_te)
        self.assertEqual(score, 2, 'wrong score')
        self.assertEqual(len(transitions_list), 4, 'wrong number of highest-scoring transitions')
        plot_array_content(query_seq, ref1, ref2, transitions_list)
    
    def test_align_2(self):
        """high ambiguity"""
        free_gap = [True, True]
        query_seq = "CTCATGCTATCTG"
        ref1      = "CTACTAACTTAGTCGA"
        ref2      = "GCTAGCTATATCGAT"
    
        score = c_align(args, query_seq, ref1, ref2, free_gap)
        transitions_list = retrieve_transitions_list(query_seq, ref1, ref2, '+', '+', 0, 0, 0, 0)
    
        for i,(ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2) in enumerate(transitions_list):
            cigar = "".join([inflate_cigar(cg_ref1), inflate_cigar(cg_gap), inflate_cigar(cg_ref2)])
            print(f"\n#{i+1}: transition {query_ts}->{query_te}, score: {score}, cigar: {cigar}")
            print_crossalignment(query_seq, ref1, ref2, cigar, fna_ref1, fa_ref2, query_ts, query_te)
        plot_array_content(query_seq, ref1, ref2, transitions_list)

    def test_align_3(self):
        """high ambiguity"""
        free_gap = [True, True]
        query_seq = "AGTAATGTGACTGGAGTTCAGACGTGTG"
        ref1      = "AGATGTGTATAAGAGACAG"
        ref2      = "GGGAAACGATGTG"

        score = c_align(args, query_seq, ref1, ref2, free_gap)
        transitions_list = retrieve_transitions_list(query_seq, ref1, ref2, '+', '+', 0, 0, 0, 0)

        for i,(ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2) in enumerate(transitions_list):
            cigar = "".join([inflate_cigar(cg_ref1), inflate_cigar(cg_gap), inflate_cigar(cg_ref2)])
            print(f"\n#{i+1}: transition {query_ts}->{query_te}, score: {score}, cigar: {cigar}")
            print_crossalignment(query_seq, ref1, ref2, cigar, fna_ref1, fa_ref2, query_ts, query_te)
        plot_array_content(query_seq, ref1, ref2, transitions_list)

if __name__ == '__main__':
    unittest.main()