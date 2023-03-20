import unittest
import os
import pickle
import pandas as pd
from crossalign import *

args = None

def write_fastq(fn, seq_id, seq, mode="w"):
    with open(fn, mode) as f:
        qual_str = "!"*len(seq)
        print(f"@{seq_id}\n{seq}\n+\n{qual_str}", file=f)

def write_fasta(fn, seq_id, seq, mode="w"):
    with open(fn, mode) as f:
        print(f">{seq_id}\n{seq}", end="", file=f)

class SimpleAlignmentTestCase(unittest.TestCase):

    def setUp(self):
        # for these test, no external files are needed
        global args
        args = parse_args(
            [
                '--reads', 'None',
                '--genome', 'None',
                '--adapter', 'None',
                '--mean', '1.15', # normally autom. determined
                '--std', '0.1'   # normally autom. determined
            ])
        init_ctypes_arrays(args)
        ct.CDLL(clib).init(args.sequence_type == 'nucl')

    def test_arguments_passed(self):
        self.assertEqual((args.match, args.mismatch, args.gap_open, args.gap_extension), (1, -3, -1, -1),
                         'arguments not parsed correctly')

    def test_align_1(self):
        free_gap = [True, True]
        query_seq = "ATGC"
        ref1      = "ATTT"
        ref2      = "GGGC"

        score = c_align(args, query_seq, ref1, ref2, free_gap)
        transitions_list = retrieve_transitions_list(query_seq, ref1, ref2, '+', '+', 0, 0, 0, 0)

        for i,(ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2) in enumerate(transitions_list):
            cigar = "".join([inflate_cigar(cg_ref1), inflate_cigar(cg_gap), inflate_cigar(cg_ref2)])
            print(f"\n#{i+1}: transition {query_ts}->{query_te}, score: {score}, cigar: {cigar}")
            print_crossalignment(query_seq, ref1, ref2, cigar, fna_ref1, fa_ref2, query_ts, query_te)
        self.assertEqual(score, 4, 'wrong score')
        self.assertEqual(len(transitions_list), 1, 'wrong number of highest-scoring transitions')
        #self.assertEqual(transitions_list[0][-1], "4=", 'wrong cigar string')
        self.assertEqual((cg_ref1, cg_gap, cg_ref2), ("2=", "", "2="))

    def test_align_2(self):
        """multiple equally high scoring alignments"""
        free_gap = [True, True]
        query_seq = "ATTC"
        ref1      = "ATTT"
        ref2      = "GCTC"

        score = c_align(args, query_seq, ref1, ref2, free_gap)
        transitions_list = retrieve_transitions_list(query_seq, ref1, ref2, '+', '+', 0, 0, 0, 0)

        for i,(ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2) in enumerate(transitions_list):
            cigar = "".join([inflate_cigar(cg_ref1), inflate_cigar(cg_gap), inflate_cigar(cg_ref2)])
            print(f"\n#{i+1}: transition {query_ts}->{query_te}, score: {score}, cigar: {cigar}")
            print_crossalignment(query_seq, ref1, ref2, cigar, fna_ref1, fa_ref2, query_ts, query_te)
        self.assertEqual(score, 4, 'wrong score')
        self.assertEqual(len(transitions_list), 2, 'wrong number of highest-scoring transitions')
        self.assertEqual(transitions_list[0][-3:], ("2=","","2="), 'wrong cigar string')
        self.assertEqual(transitions_list[1][-3:], ("3=","","1="), 'wrong cigar string')

    def test_align_3(self):
        """no "free end gap" for alignment against first reference sequence"""
        free_gap = [False, True]
        query_seq = "GACC"
        ref1      = "GATTC"
        ref2      = "ACC" 

        score = c_align(args, query_seq, ref1, ref2, free_gap)
        transitions_list = retrieve_transitions_list(query_seq, ref1, ref2, '+', '+', 0, 0, 0, 0)

        for i,(ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2) in enumerate(transitions_list):
            cigar = "".join([inflate_cigar(cg_ref1), inflate_cigar(cg_gap), inflate_cigar(cg_ref2)])
            print(f"\n#{i+1}: transition {query_ts}->{query_te}, score: {score}, cigar: {cigar}")
            print_crossalignment(query_seq, ref1, ref2, cigar, fna_ref1, fa_ref2, query_ts, query_te)
        self.assertEqual(score, 1, 'wrong score')
        self.assertEqual(len(transitions_list), 1, 'wrong number of highest-scoring transitions')
        self.assertEqual(transitions_list[0][-3:], ("2=2D1=","","1="), 'wrong cigar string')

    def test_align_4(self):
        """no "free end gap" for alignment against second reference sequence"""
        free_gap = [True, False]
        query_seq = "GACC"
        ref1      = "GACTTC"
        ref2      = "ATCC" 

        score = c_align(args, query_seq, ref1, ref2, free_gap)
        transitions_list = retrieve_transitions_list(query_seq, ref1, ref2, '+', '+', 0, 0, 0, 0)

        for i,(ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2) in enumerate(transitions_list):
            cigar = "".join([inflate_cigar(cg_ref1), inflate_cigar(cg_gap), inflate_cigar(cg_ref2)])
            print(f"\n#{i+1}: transition {query_ts}->{query_te}, score: {score}, cigar: {cigar}")
            print_crossalignment(query_seq, ref1, ref2, cigar, fna_ref1, fa_ref2, query_ts, query_te)
        self.assertEqual(score, 2, 'wrong score')
        self.assertEqual(len(transitions_list), 1, 'wrong number of highest-scoring transitions')
        self.assertEqual(transitions_list[0][-3:], ("1=","","1=1D2="), 'wrong cigar string')

    def test_align_5(self):
        """no "free end gap" for both alignments"""
        free_gap = [False, False]
        query_seq = "GACC"
        ref1      = "GAT"
        ref2      = "TCC" 

        score = c_align(args, query_seq, ref1, ref2, free_gap)
        transitions_list = retrieve_transitions_list(query_seq, ref1, ref2, '+', '+', 0, 0, 0, 0)

        for i,(ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2) in enumerate(transitions_list):
            cigar = "".join([inflate_cigar(cg_ref1), inflate_cigar(cg_gap), inflate_cigar(cg_ref2)])
            print(f"\n#{i+1}: transition {query_ts}->{query_te}, score: {score}, cigar: {cigar}")
            print_crossalignment(query_seq, ref1, ref2, cigar, fna_ref1, fa_ref2, query_ts, query_te)
        self.assertEqual(score, 1, 'wrong score')
        self.assertEqual(len(transitions_list), 1, 'wrong number of highest-scoring transitions')
        self.assertEqual(transitions_list[0][-3:], ("2=1D","","1D2="), 'wrong cigar string')

    def test_align_6(self):
        """query shorter than references (free gaps only in one direction!)"""
        free_gap = [True, True]
        query_seq = "GACC"
        ref1      = "TGA"
        ref2      = "CCT"

        score = c_align(args, query_seq, ref1, ref2, free_gap)
        transitions_list = retrieve_transitions_list(query_seq, ref1, ref2, '+', '+', 0, 0, 0, 0)

        for i,(ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2) in enumerate(transitions_list):
            cigar = "".join([inflate_cigar(cg_ref1), inflate_cigar(cg_gap), inflate_cigar(cg_ref2)])
            print(f"\n#{i+1}: transition {query_ts}->{query_te}, score: {score}, cigar: {cigar}")
            print_crossalignment(query_seq, ref1, ref2, cigar, fna_ref1, fa_ref2, query_ts, query_te)
        self.assertEqual(score, 0, 'wrong score')
        self.assertEqual(len(transitions_list), 1, 'wrong number of highest-scoring transitions')
        self.assertEqual(transitions_list[0][-3:], ("1D2=","","2=1D"), 'wrong cigar string')

    def test_align_7(self):
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

class FastqToTransitionsTestCase(unittest.TestCase):
    """Integration tests"""
    
    def setUp(self):
        self.outdir = "tests_outputs"
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def test_fastq_align_1(self):
        free_gap = [True, True]
        ref1      = "ACTGCACTACTCTAAGGACAATTACGGAGTGGACAGCCTAGCGAGAGCACCCTTCAAATGATGATAAGCACCGGCAAGCATTATTGATCAACGCAAGGATgtaatgctaaaagtccatagcacatacatcccaacctggtatgcgtacaa"
        query_seq = "ACTGCACTACTCTAAGGACAATTACGGAGTGGACAGCCTAGCGAGAGCACCCTTCAAATGATGATAAGCACCGGCAAGCATTATTGATCAACGCAAGGATCGGTGATATTAACAAAGATTCGGCACATTACTCTTGTTGGTGTGGAATCG"
        ref2      = "tatcccatcaacgtgtactcgaatattgtatattgttctcacatgaacaatttgacgatcgtttcacgctaagatgctggccacgtattaaattaatgcgCGGTGATATTAACAAAGATTCGGCACATTACTCTTGTTGGTGTGGAATCG"
        write_fastq(os.path.join(self.outdir, "reads.fastq"), "query", query_seq)
        write_fasta(os.path.join(self.outdir, "genome.fasta"), "ref1", ref1)
        write_fasta(os.path.join(self.outdir, "adapter.fasta"), "ref2", ref2)
        prefix = os.path.abspath(os.path.join(self.outdir, "1"))


        global args
        args = parse_args(
            [
                '--reads', os.path.join(self.outdir, "reads.fastq"),
                '--genome', os.path.join(self.outdir, "genome.fasta"),
                '--adapter', os.path.join(self.outdir, "adapter.fasta"),
                '--mean', '1.15', # normally autom. determined
                '--std', '0.1',   # normally autom. determined
                '--prefix', prefix,
                '--quiet'
            ])
        main(args)

        df = pd.read_pickle(prefix + '.alignment.df.pkl')
        print()
        verbose(df.explode('transitions'))
        self.assertEqual(len(df.explode('transitions')), 1)
        row = df.explode('transitions').iloc[0]
        ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2 = row.transitions
        self.assertEqual((row.subj_ref1, row.strand_ref1, int(ts), row.subj_ref2, row.strand_ref2, int(te), row.score),
                         ('ref1', '+', 100, 'ref2', '+', 100, 10), 'unexpected alignment result')

    def test_fastq_align_2(self):
        free_gap = [True, True]
        ref1      = "ACTGCACTACTCTAAGGACAATTACGGAGTGGACAGCCTAGCGAGAGCACCCTTCAAATGATGATAAGCACCGGCAAGCATTATTGATCAACGCAAGGATtaatgctaaaagtccatagcacatacatcccaacctggtatgcgtacaaa"
        query_seq = "ACTGCACTACTCTAAGGACAATTACGGAGTGGACAGCCTAGCGAGAGCACCCTTCAAATGATGATAAGCACCGGCAAGCATTATTGATCAACGCAAGGATCGGTGATATTAACAAAGATTCGGCACATTACTCTTGTTGGTGTGGAATCG"
        ref2      = "ctatcccatcaacgtgtactcgaatattgtatattgttctcacatgaacaatttgacgatcgtttcacgctaagatgctggccacgtattaaattaatgcCGGTGATATTAACAAAGATTCGGCACATTACTCTTGTTGGTGTGGAATCG"
        ref1 = ref1.translate(compl)[::-1]
        write_fastq(os.path.join(self.outdir, "reads.fastq"), "query", query_seq)
        write_fasta(os.path.join(self.outdir, "genome.fasta"), "ref1", ref1)
        write_fasta(os.path.join(self.outdir, "adapter.fasta"), "ref2", ref2)
        prefix = os.path.abspath(os.path.join(self.outdir, "2"))


        global args
        args = parse_args(
            [
                '--reads', os.path.join(self.outdir, "reads.fastq"),
                '--genome', os.path.join(self.outdir, "genome.fasta"),
                '--adapter', os.path.join(self.outdir, "adapter.fasta"),
                '--mean', '1.15', # normally autom. determined
                '--std', '0.1',   # normally autom. determined
                '--prefix', prefix,
                '--quiet'
            ])
        main(args)

        df = pd.read_pickle(prefix + '.alignment.df.pkl')
        print()
        verbose(df.explode('transitions'))
        self.assertEqual(len(df.explode('transitions')), 1)
        row = df.explode('transitions').iloc[0]
        ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2 = row.transitions
        self.assertEqual((row.subj_ref1, row.strand_ref1, int(ts), row.subj_ref2, row.strand_ref2, int(te), row.score),
                         ('ref1', '-', 50, 'ref2', '+', 100, 10), 'unexpected alignment result')

    def test_fastq_align_3(self):
        free_gap = [True, True]
        ref1      = "ACTGCACTACTCTAAGGACAATTACGGAGTGGACAGCCTAGCGAGAGCACCCTTCAAATGATGATAAGCACCGGCAAGCATTATTGATCAACGCAAGGATtaatgctaaaagtccatagcacatacatcccaacctggtatgcgtacaaa"
        query_seq = "ACTGCACTACTCTAAGGACAATTACGGAGTGGACAGCCTAGCGAGAGCACCCTTCAAATGATGATAAGCACCGGCAAGCATTATTGATCAACGCAAGGATCGGTGATATTAACAAAGATTCGGCACATTACTCTTGTTGGTGTGGAATCG"
        ref2      = "ctatcccatcaacgtgtactcgaatattgtatattgttctcacatgaacaatttgacgatcgtttcacgctaagatgctggccacgtattaaattaatgcCGGTGATATTAACAAAGATTCGGCACATTACTCTTGTTGGTGTGGAATCG"
        ref2 = ref2.translate(compl)[::-1]
        write_fastq(os.path.join(self.outdir, "reads.fastq"), "query", query_seq)
        write_fasta(os.path.join(self.outdir, "genome.fasta"), "ref1", ref1)
        write_fasta(os.path.join(self.outdir, "adapter.fasta"), "ref2", ref2)
        prefix = os.path.abspath(os.path.join(self.outdir, "3"))


        global args
        args = parse_args(
            [
                '--reads', os.path.join(self.outdir, "reads.fastq"),
                '--genome', os.path.join(self.outdir, "genome.fasta"),
                '--adapter', os.path.join(self.outdir, "adapter.fasta"),
                '--mean', '1.15', # normally autom. determined
                '--std', '0.1',   # normally autom. determined
                '--prefix', prefix,
                '--quiet'
            ])
        main(args)

        df = pd.read_pickle(prefix + '.alignment.df.pkl')
        print()
        verbose(df.explode('transitions'))
        self.assertEqual(len(df.explode('transitions')), 1)
        row = df.explode('transitions').iloc[0]
        ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2 = row.transitions
        self.assertEqual((row.subj_ref1, row.strand_ref1, int(ts), row.subj_ref2, row.strand_ref2, int(te), row.score),
                         ('ref1', '+', 100, 'ref2', '-', 50, 10), 'unexpected alignment result')

    def test_fastq_align_4(self):
        free_gap = [True, True]
        ref1      = "ACTGCACTACTCTAAGGACAATTACGGAGTGGACAGCCTAGCGAGAGCACCCTTCAAATGATGATAAGCACCGGCAAGCATTATTGATCAACGCAAGGATtaatgctaaaagtccatagcacatacatcccaacctggtatgcgtacaaa"
        query_seq = "ACTGCACTACTCTAAGGACAATTACGGAGTGGACAGCCTAGCGAGAGCACCCTTCAAATGATGATAAGCACCGGCAAGCATTATTGATCAACGCAAGGATCGGTGATATTAACAAAGATTCGGCACATTACTCTTGTTGGTGTGGAATCG"
        ref2      = "ctatcccatcaacgtgtactcgaatattgtatattgttctcacatgaacaatttgacgatcgtttcacgctaagatgctggccacgtattaaattaatgcCGGTGATATTAACAAAGATTCGGCACATTACTCTTGTTGGTGTGGAATCG"
        query_seq = query_seq.translate(compl)[::-1]
        write_fastq(os.path.join(self.outdir, "reads.fastq"), "query", query_seq)
        write_fasta(os.path.join(self.outdir, "genome.fasta"), "ref1", ref1)
        write_fasta(os.path.join(self.outdir, "adapter.fasta"), "ref2", ref2)
        prefix = os.path.abspath(os.path.join(self.outdir, "4"))


        global args
        args = parse_args(
            [
                '--reads', os.path.join(self.outdir, "reads.fastq"),
                '--genome', os.path.join(self.outdir, "genome.fasta"),
                '--adapter', os.path.join(self.outdir, "adapter.fasta"),
                '--mean', '1.15', # normally autom. determined
                '--std', '0.1',   # normally autom. determined
                '--prefix', prefix,
                '--quiet'
            ])
        main(args)

        df = pd.read_pickle(prefix + '.alignment.df.pkl')
        print()
        verbose(df.explode('transitions'))
        self.assertEqual(len(df.explode('transitions')), 1)
        row = df.explode('transitions').iloc[0]
        ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2 = row.transitions
        self.assertEqual((row.subj_ref1, row.strand_ref1, int(ts), row.subj_ref2, row.strand_ref2, int(te), row.score),
                         ('ref2', '-', 100, 'ref1', '-', 100, 10), 'unexpected alignment result')

    def test_fastq_align_5(self):
        free_gap = [True, True]
        ref1      = "ACTGCACTACTCTAAGGACAATTACGGAGTGGACAGCCTAGCGAGAGCACCCTTCAAATGATGATAAGCACCGGCAAGCATTATTGATCAA CGCAAGGATtaatgctaaaagtccatagcacatacatcccaacctggtatgcgtacaaa".replace(" ", "")
        query_seq = "ACTGCACTACTCTAAGGACAATTACGGAGTGGACAGCCTAGCGAGAGCACCCTTCAAATGATGATAAGCACCGGCAAGCATTATTGATCAAaCGC AGtATCGGTGATATTAACAAAGATTCGGCACATTACTCTTGTTGGTGTGGAATCG".replace(" ", "")
        ref2      =  "ctatcccatcaacgtgtactcgaatattgtatattgttctcacatgaacaatttgacgatcgtttcacgctaagatgctggccacgtattaaattaatgcCGGTGATATTAACAAAGATTCGGCACATTACTCTTGTTGGTGTGGAATCG"
        write_fastq(os.path.join(self.outdir, "reads.fastq"), "query", query_seq)
        write_fasta(os.path.join(self.outdir, "genome.fasta"), "ref1", ref1)
        write_fasta(os.path.join(self.outdir, "adapter.fasta"), "ref2", ref2)
        prefix = os.path.abspath(os.path.join(self.outdir, "5"))


        global args
        args = parse_args(
            [
                '--reads', os.path.join(self.outdir, "reads.fastq"),
                '--genome', os.path.join(self.outdir, "genome.fasta"),
                '--adapter', os.path.join(self.outdir, "adapter.fasta"),
                '--mean', '1.15', # normally autom. determined
                '--std', '0.1',   # normally autom. determined
                '--prefix', prefix,
                '--quiet'
            ])
        main(args)

        df = pd.read_pickle(prefix + '.alignment.df.pkl')
        print()
        verbose(df.explode('transitions'))
        self.assertEqual(len(df.explode('transitions')), 1)
        row = df.explode('transitions').iloc[0]
        ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2 = row.transitions
        self.assertEqual((row.subj_ref1, row.strand_ref1, int(ts), row.subj_ref2, row.strand_ref2, int(te)),
                         ('ref1', '+', 100, 'ref2', '+', 100))

    def test_fastq_align_6(self):
        """use blastn"""
        free_gap = [True, True]
        ref1      = "ACTGCACTACTCTAAGGACAATTACGGAGTGGACAGCCTAGCGAGAGCACCCTTCAAATGATGATAAGCACCGGCAAGCATTATTGATCAA CGCAAGGATtaatgctaaaagtccatagcacatacatcccaacctggtatgcgtacaaa".replace(" ", "")
        query_seq = "ACTGCACTACTCTAAGGACAATTACGGAGTGGACAGCCTAGCGAGAGCACCCTTCAAATGATGATAAGCACCGGCAAGCATTATTGATCAAaCGC AGtATCGGTGATATTAACAAAGATTCGGCACATTACTCTTGTTGGTGTGGAATCG".replace(" ", "")
        ref2      =  "ctatcccatcaacgtgtactcgaatattgtatattgttctcacatgaacaatttgacgatcgtttcacgctaagatgctggccacgtattaaattaatgcCGGTGATATTAACAAAGATTCGGCACATTACTCTTGTTGGTGTGGAATCG"
        write_fastq(os.path.join(self.outdir, "reads.fastq"), "query", query_seq)
        write_fasta(os.path.join(self.outdir, "genome.fasta"), "ref1", ref1)
        write_fasta(os.path.join(self.outdir, "adapter.fasta"), "ref2", ref2)
        prefix = os.path.abspath(os.path.join(self.outdir, "6"))


        global args
        args = parse_args(
            [
                '--reads', os.path.join(self.outdir, "reads.fastq"),
                '--genome', os.path.join(self.outdir, "genome.fasta"),
                '--adapter', os.path.join(self.outdir, "adapter.fasta"),
                '--adapter_alignment_tool', 'blastn',
                '--genome_alignment_tool', 'blastn',
                '--mean', '1.15', # normally autom. determined
                '--std', '0.1',   # normally autom. determined
                '--prefix', prefix,
                '--quiet'
            ])
        main(args)

        df = pd.read_pickle(prefix + '.alignment.df.pkl')
        print()
        verbose(df.explode('transitions'))
        self.assertEqual(len(df.explode('transitions')), 1)
        row = df.explode('transitions').iloc[0]
        ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2 = row.transitions
        self.assertEqual((row.subj_ref1, row.strand_ref1, int(ts), row.subj_ref2, row.strand_ref2, int(te)),
                         ('ref1', '+', 100, 'ref2', '+', 100))

    def test_fastq_align_7(self):
        """use blastn on short seqs"""
        free_gap = [True, True]
        ref1      = "AAGCATTATTGATCAA CGCAAGGATtaatgctaaaagtccatagcacata".replace(" ", "")
        query_seq = "AAGCATTATTGATCAAaCGC AGtATCGGTGATATTAACAAAGATTCGGCA".replace(" ", "")
        ref2      =  "tgctggccacgtattaaattaatgcCGGTGATATTAACAAAGATTCGGCA"
        write_fastq(os.path.join(self.outdir, "reads.fastq"), "query", query_seq)
        write_fasta(os.path.join(self.outdir, "genome.fasta"), "ref1", ref1)
        write_fasta(os.path.join(self.outdir, "adapter.fasta"), "ref2", ref2)
        prefix = os.path.abspath(os.path.join(self.outdir, "6"))


        global args
        args = parse_args(
            [
                '--reads', os.path.join(self.outdir, "reads.fastq"),
                '--genome', os.path.join(self.outdir, "genome.fasta"),
                '--adapter', os.path.join(self.outdir, "adapter.fasta"),
                '--adapter_alignment_tool', 'blastn',
                '--genome_alignment_tool', 'blastn',
                '--mean', '1.15', # normally autom. determined
                '--std', '0.1',   # normally autom. determined
                '--prefix', prefix,
                '--quiet'
            ])
        main(args)

        df = pd.read_pickle(prefix + '.alignment.df.pkl')
        print()
        verbose(df.explode('transitions'))
        self.assertEqual(len(df.explode('transitions')), 1)
        row = df.explode('transitions').iloc[0]
        ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2 = row.transitions
        self.assertEqual((row.subj_ref1, row.strand_ref1, int(ts), row.subj_ref2, row.strand_ref2, int(te)),
                         ('ref1', '+', 25, 'ref2', '+', 25))


    def test_fastq_align_7(self):
        """ambiguity due to overlap"""
        free_gap = [True, True]
        ref1      = "AAGCATTATTGATCAACGCAAGGATCaatgctaaaagtccatagcacata".replace(" ", "")
        query_seq = "AAGCATTATTGATCAACGCAAGGATCGGTGATATTAACAAAGATTCGGCA".replace(" ", "")
        ref2      = "atgctggccacgtattaaattaatTCGGTGATATTAACAAAGATTCGGCA"
        write_fastq(os.path.join(self.outdir, "reads.fastq"), "query", query_seq)
        write_fasta(os.path.join(self.outdir, "genome.fasta"), "ref1", ref1)
        write_fasta(os.path.join(self.outdir, "adapter.fasta"), "ref2", ref2)
        prefix = os.path.abspath(os.path.join(self.outdir, "7"))


        global args
        args = parse_args(
            [
                '--reads', os.path.join(self.outdir, "reads.fastq"),
                '--genome', os.path.join(self.outdir, "genome.fasta"),
                '--adapter', os.path.join(self.outdir, "adapter.fasta"),
                '--adapter_alignment_tool', 'blastn',
                '--genome_alignment_tool', 'blastn',
                '--mean', '1.15', # normally autom. determined
                '--std', '0.1',   # normally autom. determined
                '--prefix', prefix,
                '--quiet'
            ])
        main(args)

        df = pd.read_pickle(prefix + '.alignment.df.pkl')
        print()
        verbose(df.explode('transitions'))
        self.assertEqual(len(df.explode('transitions')), 3)
        #row = df.explode('transitions').iloc[0]
        #ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2 = row.transitions
        #self.assertEqual((row.subj_ref1, row.strand_ref1, int(ts), row.subj_ref2, row.strand_ref2, int(te)),
        #                 ('ref1', '+', 25, 'ref2', '+', 25))

    def test_fastq_align_8(self):
        """resolve ambiguity by using site of interest"""
        free_gap = [True, True]
        ref1      = "AAGCATTATTGATCAACGCAAGGATCaatgctaaaagtccatagcacata".replace(" ", "")
        query_seq = "AAGCATTATTGATCAACGCAAGGATCGGTGATATTAACAAAGATTCGGCA".replace(" ", "")
        ref2      = "atgctggccacgtattaaattaatTCGGTGATATTAACAAAGATTCGGCA"
        write_fastq(os.path.join(self.outdir, "reads.fastq"), "query", query_seq)
        write_fasta(os.path.join(self.outdir, "genome.fasta"), "ref1", ref1)
        write_fasta(os.path.join(self.outdir, "adapter.fasta"), "ref2", ref2)
        prefix = os.path.abspath(os.path.join(self.outdir, "8"))
        with open(os.path.join(self.outdir, "soi.tsv"), "w") as f:
            print("ref1\t+\t25", file=f)


        global args
        args = parse_args(
            [
                '--reads', os.path.join(self.outdir, "reads.fastq"),
                '--genome', os.path.join(self.outdir, "genome.fasta"),
                '--adapter', os.path.join(self.outdir, "adapter.fasta"),
                '--sites_of_interest', os.path.join(self.outdir, "soi.tsv"),
                '--adapter_alignment_tool', 'blastn',
                '--genome_alignment_tool', 'blastn',
                '--mean', '1.15', # normally autom. determined
                '--std', '0.1',   # normally autom. determined
                '--prefix', prefix,
                '--quiet'
            ])
        main(args)

        df = pd.read_pickle(prefix + '.alignment.df.pkl')
        print()
        verbose(df.explode('transitions'))
        self.assertEqual(len(df.explode('transitions')), 1)
        row = df.explode('transitions').iloc[0]
        ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2 = row.transitions
        self.assertEqual((row.subj_ref1, row.strand_ref1, int(ts), row.subj_ref2, row.strand_ref2, int(te)),
                         ('ref1', '+', 25, 'ref2', '+', 25))

    def test_fastq_align_9(self):
        """all possible combinations of transition order and orientations"""
        free_gap = [True, True]
        ref1        = "AAGCATTATTGATCAACGCAAGGAT gaatgctaaaagtccatagcacata".replace(" ", "")
        ref2        = "atgctggccacgtattaaattaatg CGGTGATATTAACAAAGATTCGGCA".replace(" ", "")
        
        query_12_pp = "AAGCATTATTGATCAACGCAAGGAT CGGTGATATTAACAAAGATTCGGCA".replace(" ", "")
        query_12_pm = "AAGCATTATTGATCAACGCAAGGAT cattaatttaatacgtggccagcat".replace(" ", "")
        query_12_mp = "tatgtgctatggacttttagcattc CGGTGATATTAACAAAGATTCGGCA".replace(" ", "")
        query_12_mm = "tatgtgctatggacttttagcattc cattaatttaatacgtggccagcat".replace(" ", "")
        query_21_pp = "atgctggccacgtattaaattaatg gaatgctaaaagtccatagcacata".replace(" ", "")
        query_21_pm = "atgctggccacgtattaaattaatg ATCCTTGCGTTGATCAATAATGCTT".replace(" ", "")
        query_21_mp = "TGCCGAATCTTTGTTAATATCACCG gaatgctaaaagtccatagcacata".replace(" ", "")
        query_21_mm = "TGCCGAATCTTTGTTAATATCACCG ATCCTTGCGTTGATCAATAATGCTT".replace(" ", "")

        write_fastq(os.path.join(self.outdir, "reads.fastq"), "12_pp", query_12_pp)
        write_fastq(os.path.join(self.outdir, "reads.fastq"), "12_pm", query_12_pm, 'a')
        write_fastq(os.path.join(self.outdir, "reads.fastq"), "12_mp", query_12_mp, 'a')
        write_fastq(os.path.join(self.outdir, "reads.fastq"), "12_mm", query_12_mm, 'a')
        write_fastq(os.path.join(self.outdir, "reads.fastq"), "21_pp", query_21_pp, 'a')
        write_fastq(os.path.join(self.outdir, "reads.fastq"), "21_pm", query_21_pm, 'a')
        write_fastq(os.path.join(self.outdir, "reads.fastq"), "21_mp", query_21_mp, 'a')
        write_fastq(os.path.join(self.outdir, "reads.fastq"), "21_mm", query_21_mm, 'a')
        write_fasta(os.path.join(self.outdir, "genome.fasta"), "ref1", ref1)
        write_fasta(os.path.join(self.outdir, "adapter.fasta"), "ref2", ref2)
        prefix = os.path.abspath(os.path.join(self.outdir, "9"))

        global args
        args = parse_args(
            [
                '--reads', os.path.join(self.outdir, "reads.fastq"),
                '--genome', os.path.join(self.outdir, "genome.fasta"),
                '--adapter', os.path.join(self.outdir, "adapter.fasta"),
                '--adapter_alignment_tool', 'blastn',
                '--genome_alignment_tool', 'blastn',
                '--mean', '1.15', # normally autom. determined
                '--std', '0.1',   # normally autom. determined
                '--prefix', prefix,
                '--quiet'
            ])
        main(args)

        df = pd.read_pickle(prefix + '.alignment.df.pkl')
        print()
        verbose(df.explode('transitions'))
        self.assertEqual(len(df.explode('transitions')), 8)
        for i,row in df.explode('transitions').iterrows():
            ts, te, fna_ref1, fa_ref2, query_ts, query_te, cg_ref1, cg_gap, cg_ref2 = row.transitions
            self.assertEqual((int(ts), int(te)),
                             (25, 25))


if __name__ == '__main__':
    unittest.main()