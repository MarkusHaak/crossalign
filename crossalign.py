import logging
import sys
import os, argparse, shutil
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from scipy.stats import norm
import ctypes as ct
import re
from pandarallel import pandarallel
from itertools import groupby, count
from collections import deque
from tqdm import tqdm
import pysam

# translate for building complement of a DNA sequence
compl = str.maketrans('ATGCNSWRYKMatgcnswrykm', 'TACGNWSYRMKtacgnwsyrmk')
cg_pat = re.compile("(\d+)([=XID])")
plt.rcParams.update({'font.size': 12})

# initialize the necessary data structures
clib = os.path.join(os.path.dirname(os.path.realpath(__file__)), "align.so")

nd_pp = np.ctypeslib.ndpointer(dtype=np.uintp, ndim=1, flags='C_CONTIGUOUS')
align = ct.CDLL(clib).align
align.argtypes = [ct.c_char_p, ct.c_char_p, 
                  ct.c_short, ct.c_short, ct.c_short,
                  ct.c_short, ct.c_short, ct.c_short, ct.c_short,
                  ct.c_bool, ct.c_bool,
                  nd_pp, nd_pp, nd_pp, nd_pp, nd_pp, nd_pp, nd_pp]
align.restype = ct.c_int

backtrace = ct.CDLL(clib).backtrace
backtrace.argtypes = [nd_pp, nd_pp, nd_pp, nd_pp, nd_pp, nd_pp,
                      ct.c_short, ct.c_short, ct.c_short,
                      nd_pp, nd_pp, nd_pp]
backtrace.restype = ct.c_int

get_cigar = ct.CDLL(clib).get_cigar
get_cigar.argtypes = [nd_pp, nd_pp, nd_pp, nd_pp, nd_pp, nd_pp,
                      ct.c_short, ct.c_short, ct.c_short,
                      nd_pp, ct.c_short, ct.c_short, ct.c_short, ct.c_short, ct.c_char_p]
get_cigar.restype = ct.c_int

get_transitions = ct.CDLL(clib).get_transitions
get_transitions.argtypes = [nd_pp, ct.c_short, ct.c_short, ct.c_short,
                            nd_pp, nd_pp, nd_pp, nd_pp, nd_pp]
get_transitions.restype = ct.c_int

class ArgHelpFormatter(argparse.HelpFormatter):
    '''
    Formatter adding default values to help texts.
    '''
    def __init__(self, prog):
        super().__init__(prog)

    ## https://stackoverflow.com/questions/3853722
    #def _split_lines(self, text, width):
    #   if text.startswith('R|'):
    #       return text[2:].splitlines()  
    #   # this is the RawTextHelpFormatter._split_lines
    #   return argparse.HelpFormatter._split_lines(self, text, width)

    def _get_help_string(self, action):
        text = action.help
        if  action.default is not None and \
            action.default != argparse.SUPPRESS and \
            'default:' not in text.lower():
            text += ' (default: {})'.format(action.default)
        return text

def get_args(args=None):
    parser = argparse.ArgumentParser(description='Estimates read starts (transposase insertion sites) for ONT rapid libraries',
                                     formatter_class=ArgHelpFormatter, 
                                     add_help=False)

    main_group = parser.add_argument_group('Main Arguments')
    main_group.add_argument('--reads',
                            required=True,
                            nargs='+',
                            help='fastq files or path to directories containing fastq files (recursion depth 1)')
    main_group.add_argument('--genome',
                            required=True,
                            help='Fasta file containing the genomic sequences that is searched for insertion sites.')
    main_group.add_argument('--adapter',
                            required=True,
                            help='Transposon Y adapter sequence')
    main_group.add_argument('--sites_of_interest',
                            help='''Force all alignments of potential transitions to align exactly to the given sites, leaving 
                            only the transition location of the second reference as a free variable. The required format of this
                            file is tab-separated values of subject names and 0-based sites as first and second column values
                            , respectively, without header.''')
    main_group.add_argument('--genome_alignment',
                            help='Alignment file (.paf with cigar string or sorted and indexed .bam) containing the reads vs. genome reference mapping results.')
    main_group.add_argument('--adapter_alignment',
                            help='Alignment file (.paf with cigar string or sorted and indexed .bam) containing the reads vs. adapter reference mapping results.')
    main_group.add_argument('--prefix',
                            help="filename prefix (optionally with absolute path) for generated output files.",
                            default="crossalign")
    main_group.add_argument('--plot',
                            help='plot results of gaussian sequencing error approximation',
                            action='store_true')
    main_group.add_argument('--circular',
                            action="store_true")
    main_group.add_argument('--strip',
                            help="minimal number of bases stripped from local alignments to cope with coincidently identical sequences",
                            type=int,
                            default=5)
    main_group.add_argument('--wordsize',
                            help='''minimal number of consecutive terminal matches that need to remain after trimming bases from
                            local alignments.''',
                            type=int,
                            default=5)
    main_group.add_argument('--max_dist',
                            help='''max distance between an adapter and a genome alignment to perform pairwise-alignment 
                                 in order to identify the exact transition point''',
                            type=int,
                            default=200)
    main_group.add_argument('--sequence_type',
                            default='nucl',
                            choices=['nucl', 'aa'])
    main_group.add_argument('--processes',
                            type=int,
                            default=6)
    main_group.add_argument('--progress',
                            help='show progress bars',
                            action='store_true')
    main_group.add_argument('-q','--quiet',
                            help='limit output to warning and error messages',
                            action='store_true')
    
    align_group = parser.add_argument_group('Alignment Arguments')
    align_group.add_argument('--match',
                             help='match score',
                             type=int,
                             default=1)
    align_group.add_argument('--mismatch',
                             help='mismatch penalty',
                             type=int,
                             default=-2)
    align_group.add_argument('--gap_open',
                             help='gap open penalty',
                             type=int,
                             default=-1)
    align_group.add_argument('--gap_extension',
                             help='gap extension penalty',
                             type=int,
                             default=-1)

    filter_group = parser.add_argument_group('Filter Arguments')
    filter_group.add_argument('--mean',
                              help='mean of per-base difference of actual sequence length from read length',
                              type=float)
    filter_group.add_argument('--std',
                              help='standard deviation of per-base difference of actual sequence length from read length',
                              type=float)
    filter_group.add_argument('--min_adapter_blen',
                              help="min produced alignment length (including errors)",
                              type=int,
                              default=1)
    filter_group.add_argument('--min_genome_blen',
                              help="min produced alignment length (including errors)",
                              type=int,
                              default=1)

    help_group = parser.add_argument_group('Help')
    help_group.add_argument('-h', '--help', 
                            action='help', 
                            default=argparse.SUPPRESS,
                            help='Show this help message and exit.')
    if args:
        return parser.parse_args(args)
    return parser.parse_args()

def run_minimap2(ref_fn, fq_fn, paf_fn):
    cmd ='minimap2 -x map-ont -c --eqx --secondary=no -t 4 {} {} >{} 2> /dev/null'.format(ref_fn, fq_fn, paf_fn)
    return os.system(cmd)

def parse_paf(fn, cigar=False):
    usecols = list(range(9))
    names = ["qid", "qlen", "qst", "qen", "strand", "subj", 
             "slen", "sst", "sen"]
    dtype = {"qid": str, "qlen": np.int32, "qst": np.int32, 
             "qen": np.int32, "strand": str, "subj": str,
             "slen": np.int32, "sst": np.int32, "sen": np.int32}
    converters = {}
    if cigar:
        usecols.append(22)
        names.append('cg')
        converters['cg'] = lambda x: x.split(':')[-1].replace('M', '=')
    return pd.read_csv(fn, sep='\t', header=None,
                       usecols=usecols,
                       names=names,
                       dtype=dtype,
                       converters=converters)

def parse_bam(fn, reads):
    f = pysam.AlignmentFile(fn, 'rb')
    names = ["qid", "qlen", "qst", "qen", "strand", "subj", 
             "slen", "sst", "sen", "cg"]
    data = []
    for seg in f.fetch():
        m = re.fullmatch(r"(\d+[HS])*((?:\d+[^HS\d])*)(\d+[HS])*", seg.cigarstring)
        qid = seg.query_name
        qlen = sum([j for i,j in seg.cigartuples if i != 2]) #len(readseg.loc[qid, 'seq'])##seg.template_length #seg.query_length
        qst = seg.query_alignment_start
        qen = seg.query_alignment_end
        if (seg.is_reverse == True and m.group(3)) or \
           (seg.is_reverse == False and m.group(1)):
            if seg.is_reverse == True and 'H' in m.group(3):
                clipped = int(m.group(3).rstrip('H'))
                qst += clipped
                qen += clipped
            elif seg.is_reverse == False and 'H' in m.group(1):
                clipped = int(m.group(1).rstrip('H'))
                qst += clipped
                qen += clipped
            elif seg.is_reverse == True and 'S' in m.group(3):
                clipped = int(m.group(3).rstrip('S'))
                qen = clipped + (qen - qst)
                qst = clipped

        strand = '-' if seg.is_reverse else '+'
        subj = seg.reference_name
        slen = f.get_reference_length(subj) #seg.reference_length # "aligned length of the read on the reference genome" --> could be wrong field
        sst = seg.reference_start
        sen = seg.reference_end
        cg = m.group(2).replace('M', '=')
        data.append( (qid, qlen, qst, qen, strand, subj, slen, sst, sen, cg) )
    return pd.DataFrame(data, columns=names)

def print_alignment(query, ref1, ref2, cigar, r1_fna, r2_fa, query_ts, query_te, width=50):
    ref1_ss = []
    ref1_op = []
    quer_ss = []
    ref2_op = []
    ref2_ss = []

    query_loc = 0
    ref1_loc = 0
    ref2_loc = r2_fa
    for k,op in enumerate(cigar):
        if op == '=':
            if query_loc < query_te or ref1_loc < r1_fna:
                ref1_ss.append(ref1[ref1_loc])
                ref1_op.append('|')
                ref1_loc += 1
            else:
                if ref1_loc < len(ref1):
                    ref1_ss.append(ref1[ref1_loc].lower())
                    ref1_loc += 1
                else:
                    ref1_ss.append(' ')
                ref1_op.append(' ')
            if query_loc >= query_ts and ref1_loc >= r1_fna:
                ref2_ss.append(ref2[ref2_loc])
                ref2_op.append('|')
                ref2_loc += 1
            else:
                #ref2_ss.append(' ')
                ref2_op.append(' ')
            quer_ss.append(query[query_loc])
            query_loc += 1
        elif op == 'X':
            if query_loc < query_te or ref1_loc < r1_fna:
                ref1_ss.append(ref1[ref1_loc])
                ref1_op.append('X')
                ref1_loc += 1
            else:
                if ref1_loc < len(ref1):
                    ref1_ss.append(ref1[ref1_loc].lower())
                    ref1_loc += 1
                else:
                    ref1_ss.append(' ')
                ref1_op.append(' ')
            if query_loc >= query_ts and ref1_loc >= r1_fna:
                ref2_ss.append(ref2[ref2_loc])
                ref2_op.append('X')
                ref2_loc += 1
            else:
                #ref2_ss.append(' ')
                ref2_op.append(' ')
            quer_ss.append(query[query_loc])
            query_loc += 1
        elif op == 'I':
            if query_loc < query_te or ref1_loc < r1_fna:
                ref1_ss.append('-')
                ref1_op.append('I')
                #ref1_loc += 1
            else:
                if ref1_loc < len(ref1):
                    ref1_ss.append(ref1[ref1_loc].lower())
                    ref1_loc += 1
                else:
                    ref1_ss.append(' ')
                ref1_op.append(' ')
            if query_loc >= query_ts and ref1_loc >= r1_fna:
                ref2_ss.append('-')
                ref2_op.append('I')
                #ref2_loc += 1
            else:
                #ref2_ss.append(' ')
                ref2_op.append(' ')
            quer_ss.append(query[query_loc])
            query_loc += 1
        elif op == 'D':
            if query_loc < query_te or ref1_loc < r1_fna:
                ref1_ss.append(ref1[ref1_loc]) if ref1_loc < len(ref1) else ref1_ss.append("#")
                ref1_op.append('D')
                ref2_op.append(' ')
                ref1_loc += 1
            else:
                if ref1_loc < len(ref1):
                    ref1_ss.append(ref1[ref1_loc].lower())
                    ref1_loc += 1
                else:
                    ref1_ss.append(' ')
                ref1_op.append(' ')
                if query_loc >= query_te:
                    ref2_ss.append(ref2[ref2_loc])
                    ref2_op.append('D')
                    ref2_loc += 1
                else:
                    #ref2_ss.append(' ')
                    ref2_op.append(' ')
            quer_ss.append('-')
            #query_loc += 

    to_prepend = []
    for i in range(1,len(ref2_op) - len(ref2_ss) +1):
        ref2_loc = r2_fa - i
        if ref2_loc >= 0:
            to_prepend.append(ref2[ref2_loc].lower())
        else:
            to_prepend.append(' ')
    ref2_ss = to_prepend[::-1] + ref2_ss

    for i in range(int(np.round(len(quer_ss) / width + .5))):
        for line_start, l in [('ref1:', ref1_ss), ('', ref1_op), ('query:', quer_ss), ('', ref2_op), ('ref2:', ref2_ss)]:
            print("{:<7}{}".format(line_start, "".join(l[i*width:(i+1)*width])))
        print()

def inflate_cigar(cigar):
    return ''.join([int(m.group(1)) * m.group(2) for m in re.finditer(r'(\d+)(\D)', cigar)])

def verbose(df):
    for i, row in df.iterrows():
        ref1, ref2 = get_reference_sequences(row)
        query_seq = reads.loc[row.rid].seq[int(row.qst) : int(row.qen)]
        ts, te, fna_ref1, fa_ref2, query_ts, query_te, cigar = row.transitions
        cigar = inflate_cigar(cigar)
        print('"{}" ({}) @{} -> "{}" ({}) @{}'.format(row.subj_ref1, row.strand_ref1, ts, row.subj_ref2, row.strand_ref2, te))
        print_alignment(query_seq, ref1, ref2, cigar, fna_ref1, fa_ref2, query_ts, query_te, width=50)

def traverse_cg(trimmed_query, trimmed_subject, bases, op):
    if op == "=" and bases >= args.wordsize and trimmed_query >= args.strip:
        return trimmed_query, trimmed_subject, True
    if op == "=" or op == 'X':
        trimmed_query += bases
        trimmed_subject += bases
    elif op == 'I':
        trimmed_query += bases
    elif op == 'D':
        trimmed_subject += bases
    else:
        logger.error('ERROR: unknown CIGAR Operation:', op)
        exit(1)
    return trimmed_query, trimmed_subject, False

def traverse_cg_forwards(cg):
    trimmed_query = 0
    trimmed_subject = 0
    for m in re.finditer(cg_pat, cg):
        trimmed_query, trimmed_subject, done = traverse_cg(trimmed_query, trimmed_subject, int(m.group(1)), m.group(2))
        if done:
            return trimmed_query, trimmed_subject
    return np.nan, np.nan

def traverse_cg_backwards(cg):
    trimmed_query = 0
    trimmed_subject = 0
    cg_list = re.findall(cg_pat, cg)
    for bases,op in reversed(cg_list):
        trimmed_query, trimmed_subject, done = traverse_cg(trimmed_query, trimmed_subject, int(bases), op)
        if done:
            return trimmed_query, trimmed_subject
    return np.nan, np.nan

def fix_qst(row):
    if row.trans_order == 1.:
        cg = row.cg_ad
        qen = row.qen_ad
        strand = row.strand_ad
    else:
        cg = row.cg_gn
        qen = row.qen_gn
        strand = row.strand_gn

    if strand == '+':
        trimmed_query, trimmed_subject = traverse_cg_backwards(cg)
        if trimmed_subject:
            return qen - trimmed_query, trimmed_subject
    else:
        trimmed_query, trimmed_subject = traverse_cg_forwards(cg)
        if trimmed_subject:
            return qen - trimmed_query, trimmed_subject
    # unsuccessful
    return np.nan, np.nan

def fix_qen(row):
    if row.trans_order == 1.:
        cg = row.cg_gn
        qst = row.qst_gn
        strand = row.strand_gn
    else:
        cg = row.cg_ad
        qst = row.qst_ad
        strand = row.strand_ad

    if strand == '+':
        trimmed_query, trimmed_subject = traverse_cg_forwards(cg)
        if trimmed_subject:
            return qst + trimmed_query, trimmed_subject
    else:
        trimmed_query, trimmed_subject = traverse_cg_backwards(cg)
        if trimmed_subject:
            return qst + trimmed_query, trimmed_subject
    # unsuccessful
    return np.nan, np.nan

def set_qst_and_qen(df, sel, fix=True):
    df['trim_ref1'], df['trim_ref2'] = 0, 0
    cg_ad = df.cg_ad.str.replace('D|I', 'X')
    cg_gn = df.cg_gn.str.replace('D|I', 'X')
    # set query start and end for trivial cases
    for trans_order, strand_ref1, strand_ref2, cg_ref1, cg_ref2, ref1_qen, ref2_qst in [( 1., 'strand_ad', 'strand_gn', cg_ad, cg_gn, 'qen_ad', 'qst_gn'),
                                                                                        (-1., 'strand_gn', 'strand_ad', cg_gn, cg_ad, 'qen_gn', 'qst_ad')]:
        logger.info("setting {} qst".format(trans_order))
        sel_ = ((df[strand_ref1] == '+') & (cg_ref1.str.rstrip('=').str.rsplit('X', n=1).str[-1].astype(np.float32) >= (args.wordsize + args.strip))) | \
               ((df[strand_ref1] == '-') & (cg_ref1.str.split('=', n=1).str[0].astype(np.float32) >= (args.wordsize + args.strip)))
        #sel_ = ((df[strand_ref1] == '+') & (df[cg_ref1].str.extract(r'(\d)=$').astype(np.float32) >= (args.wordsize + args.strip))) | \
        #       ((df[strand_ref1] == '-') & (df[cg_ref1].str.extract(r'^(\d)=').astype(np.float32) >= (args.wordsize + args.strip)))
        df.loc[sel & (df.trans_order == trans_order) & sel_ , 'qst'] = df.loc[sel & (df.trans_order == trans_order) & sel_, ref1_qen]
        logger.info("setting {} qen".format(trans_order))
        sel_ = ((df[strand_ref2] == '+') & (cg_ref2.str.split('=', n=1).str[0].astype(np.float32) >= (args.wordsize + args.strip))) | \
               ((df[strand_ref2] == '-') & (cg_ref2.str.rstrip('=').str.rsplit('X', n=1).str[-1].astype(np.float32) >= (args.wordsize + args.strip))) 
        #       ((df[strand_ref2] == '-') & (df[cg_ref2].str.rstrip('=').str.rsplit('X', n=1).str[-1].astype(np.float32) >= (args.wordsize + args.strip)))
        df.loc[sel & (df.trans_order == trans_order) & sel_ , 'qen'] = df.loc[sel & (df.trans_order == trans_order) & sel_, ref2_qst]
    
    sel_ = sel & df.qst.notnull()
    df.loc[sel_, 'qst'] -= args.strip
    df.loc[sel_, 'trim_ref1'] += args.strip
    sel_ = sel & df.qen.notnull()
    df.loc[sel_, 'qen'] += args.strip
    df.loc[sel_, 'trim_ref2'] += args.strip
    
    # for query start and end for non-trivial cases
    if fix:
        logger.info('need to fix {} query starts'.format(sum(sel & df.qst.isnull())))
        df.loc[sel & df.qst.isnull(), ['qst', 'trim_ref1']] = pd.DataFrame(
            df.loc[sel & df.qst.isnull()].parallel_apply(lambda row: fix_qst(row), axis=1).values.tolist(), 
            index=df.loc[sel & df.qst.isnull()].index, columns=['qst', 'trim_ref1']
        )
        logger.info('need to fix {} query ends'.format(sum(sel & df.qen.isnull())))
        df.loc[sel & df.qen.isnull(), ['qen', 'trim_ref2']] = pd.DataFrame(
            df.loc[sel & df.qen.isnull()].parallel_apply(lambda row: fix_qen(row), axis=1).values.tolist(), 
            index=df.loc[sel & df.qen.isnull()].index, columns=['qen', 'trim_ref2']
        )
    return df

def sequence_length_stats(df):
    dfs = []
    # add length of read seq that aligns as "qalen"
    dfs.append((df[df.subj.notnull()].qen - df[df.subj.notnull()].qst).to_frame(name='qalen'))
    # add length of genome seq that aligns as "salen"
    dfs.append((df[df.subj.notnull()].sen - df[df.subj.notnull()].sst).abs().to_frame(name='salen'))

    stats = pd.concat(dfs, axis=1).fillna(0.)
    data = (stats.salen - stats.qalen)/stats.qalen
    args.mean, args.std = norm.fit(data)
    if args.plot:
        fig, ax = plt.subplots(figsize=(8,6))
        ax.hist(data, bins=100, density=True)
        xmin, xmax = ax.get_xlim()
        x = np.linspace(xmin, xmax, 200)
        y = norm.pdf(x, args.mean, args.std)
        ax.plot(x, y)
        plt.show()

def count_iter_items(iterable):
        counter = count()
        deque(zip(iterable, counter), maxlen=0)  # (consume at C speed)
        return next(counter)

def c_align_row(row):
    if pd.isnull(row.max_ref_len):
        return row

    free_gap = [True, True]
    if args.sites_of_interest:
        # align only those reads that span a site of interest
        hit = False
        if (row.subj_ref1, row.strand_ref1) in soi.index:
            d = soi.loc[[(row.subj_ref1, row.strand_ref1)]]
            d = d.loc[(row.sst_ref1 <= d.site) & (d.site <= row.sen_ref1)]
            if len(d) == 1:
                if row.strand_ref1 == '+':
                    row.sen_ref1 = d.site[0]
                else:
                    row.sst_ref1 = d.site[0]
                hit ^= True
                free_gap[0] = False
        opp_strand_ref2 = '+' if row.strand_ref2 == '-' else '-'
        if (row.subj_ref2, opp_strand_ref2) in soi.index:
            d = soi.loc[[(row.subj_ref2, opp_strand_ref2)]]
            d = d.loc[(row.sst_ref2 <= d.site) & (d.site <= row.sen_ref2)]
            if len(d) == 1:
                if row.strand_ref2 == '+':
                    row.sst_ref2 = d.site[0]
                else:
                    row.sen_ref2 = d.site[0]
                hit ^= True
                free_gap[1] = False
        if hit == False:
            return row

    # determine query and reference seqeunces
    query_seq = reads.loc[row.rid].seq[int(row.qst) : int(row.qen)]
    ref1, ref2 = get_reference_sequences(row)

    logger.debug("{}, {}, {}, {}".format(query_seq, ref1, ref2, free_gap))
    qlen, s1len, s2len = len(query_seq), len(ref1), len(ref2)
    query, subj = query_seq.encode("utf8"), (ref1+ref2).encode("utf8")
    assert align(query, subj, qlen, s1len, s2len,
                 args.match, args.mismatch, args.gap_open, args.gap_extension,
                 *free_gap,
                 scores_pp, *ops_pp) == 0, "alignment failed"
    score = scores[s1len+s2len, qlen]

    assert backtrace(*ops_pp, qlen, s1len, s2len,
                     align_ends_pp0, align_ends_pp1, reachable_pp) == 0 , \
           "failed to trace back the alignment's end points"
    assert get_transitions(ops_pp[0], qlen, s1len, s2len,
                           align_ends_pp0, align_ends_pp1, reachable_pp, 
                           transitions_pp0, transitions_pp1) == 0, \
           "failed to determine transition sites based on alignment endpoints"

    transitions_list = []
    for fna_ref1, fa_ref2 in zip(*np.where(transitions[0,:s1len+1, :s2len+1] >= 0)):
        query_ts = transitions[0, fna_ref1, fa_ref2]
        query_te = transitions[1, fna_ref1, fa_ref2]
        assert get_cigar(*ops_pp, qlen, s1len, s2len,
                         reachable_pp, fna_ref1, fa_ref2, query_ts, query_te,
                         cigarbuffer) == 0, \
               "failed to retrieve the CIGAR string for alignment"
        cigar = cigarbuffer.value.decode('utf-8')
        cigar = "".join(["{}{}".format(count_iter_items(g), k) for k,g in groupby(cigar)])

        if row.strand_ref1 == '+':
            ts = row.sst_ref1 + fna_ref1
        else:
            ts = row.sen_ref1 - fna_ref1
        if row.strand_ref2 == '+':
            te = row.sst_ref2 + fa_ref2
        else:
            te = row.sen_ref2 - fa_ref2

        transitions_list.append( (ts, te, fna_ref1, fa_ref2, query_ts, query_te, cigar) )

    row['score'] = score
    row['transitions'] = transitions_list
    return row

def set_references(df, sel):
    sel_ = sel & (df.trans_order == 1.)
    if sel_.any():
        df.loc[sel_, 'subj_ref1']   = df.loc[sel_, 'subj_ad']
        df.loc[sel_, 'strand_ref1'] = df.loc[sel_, 'strand_ad']
        df.loc[sel_, 'sst_ref1']    = df.loc[sel_, 'sst_ad']
        df.loc[sel_, 'sen_ref1']    = df.loc[sel_, 'sen_ad']
        df.loc[sel_, 'subj_ref2']   = df.loc[sel_, 'subj_gn']
        df.loc[sel_, 'strand_ref2'] = df.loc[sel_, 'strand_gn']
        df.loc[sel_, 'sst_ref2']    = df.loc[sel_, 'sst_gn']
        df.loc[sel_, 'sen_ref2']    = df.loc[sel_, 'sen_gn']
    sel_ = sel & (df.trans_order == -1.)
    if sel_.any():
        df.loc[sel_, 'subj_ref1']   = df.loc[sel_, 'subj_gn']
        df.loc[sel_, 'strand_ref1'] = df.loc[sel_, 'strand_gn']
        df.loc[sel_, 'sst_ref1']    = df.loc[sel_, 'sst_gn']
        df.loc[sel_, 'sen_ref1']    = df.loc[sel_, 'sen_gn']
        df.loc[sel_, 'subj_ref2']   = df.loc[sel_, 'subj_ad']
        df.loc[sel_, 'strand_ref2'] = df.loc[sel_, 'strand_ad']
        df.loc[sel_, 'sst_ref2']    = df.loc[sel_, 'sst_ad']
        df.loc[sel_, 'sen_ref2']    = df.loc[sel_, 'sen_ad']

    sel_ = sel & (df.strand_ref1 == '+')
    if sel_.any():
        df.loc[sel_, 'sst_ref1']    = df[sel_].sen_ref1 - df[sel_].trim_ref1
        df.loc[sel_, 'sen_ref1']    = df[sel_].sst_ref1 + df[sel_].max_ref_len
    sel_ = sel & (df.strand_ref1 == '-')
    if sel_.any():
        df.loc[sel_, 'sen_ref1']    =  df[sel_].sst_ref1 + df[sel_].trim_ref1
        df.loc[sel_, 'sst_ref1']    = (df[sel_].sen_ref1 - df[sel_].max_ref_len).clip(lower=0)
    sel_ = sel & (df.strand_ref2 == '+')
    if sel_.any():
        df.loc[sel_, 'sen_ref2']    =  df[sel_].sst_ref2 + df[sel_].trim_ref2
        df.loc[sel_, 'sst_ref2']    = (df[sel_].sen_ref2 - df[sel_].max_ref_len).clip(lower=0)
    sel_ = sel & (df.strand_ref2 == '-')
    if sel_.any():
        df.loc[sel_, 'sst_ref2']    = df[sel_].sen_ref2 - df[sel_].trim_ref2
        df.loc[sel_, 'sen_ref2']    = df[sel_].sst_ref2 + df[sel_].max_ref_len

    return df

def get_reference_sequences(row):
    if row.trans_order == 1.:
        ref1 = adapter.loc[row.subj_ref1].seq[int(row.sst_ref1) : int(row.sen_ref1)]
        ref2 = genome.loc[row.subj_ref2].seq[int(row.sst_ref2) : int(row.sen_ref2)]
    else:
        ref1 = genome.loc[row.subj_ref1].seq[int(row.sst_ref1) : int(row.sen_ref1)]
        ref2 = adapter.loc[row.subj_ref2].seq[int(row.sst_ref2) : int(row.sen_ref2)]
    if row.strand_ref1 == '-':
        ref1 = ref1.translate(compl)[::-1]
    if row.strand_ref2 == '-':
        ref2 = ref2.translate(compl)[::-1]
    return ref1, ref2

def plot_norm_score_distribution(df, title, nbins=50):
    x, y = df.qlen, df.norm_score
    
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(figsize=(5, 5))

    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='in', top=True, right=True)
    ax_histx = plt.axes(rect_histx)
    ax_histx.tick_params(direction='in', labelbottom=False)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=False)

    ax_scatter.scatter(x, y, s=1., alpha=0.05)

    ax_histx.hist(x, bins=nbins)
    ax_histy.hist(y, bins=nbins, orientation='horizontal')

    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())
    
    ax_scatter.set_xlabel("distance between alignments / nt")
    ax_scatter.set_ylabel("qlen normalized alignment score")
    ax_histx.set_ylabel("bin count")
    ax_histy.set_xlabel("bin count")
    
    ax_histx.set_title(title)
    
    plt.show()


if __name__ == '__main__':
    args = get_args()

    if args.quiet:
        logging.basicConfig(stream=sys.stdout, level=logging.WARNING, 
                            format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S')
    else:
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                            format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S')
    logger = logging.getLogger('main')

    pandarallel.initialize(nb_workers=args.processes, progress_bar=(args.progress and not args.quiet))

    # read sequence data
    fq_files = []
    for entry in args.reads:
        if os.path.isfile(entry) and (entry.endswith(".fastq") or entry.endswith(".fq")):
            fq_files.append(entry)
        else:
            fq_files.extend([os.path.join(entry, f) for f in os.listdir(entry) if os.path.isfile(os.path.join(entry, f)) \
                and (f.endswith(".fastq") or f.endswith(".fq"))])

    if args.sites_of_interest:
        logger.info(" - reject potential transitions that do not span a site of interest")
        soi = pd.read_csv(args.sites_of_interest, header=None, sep='\t', names=['subj', 'strand', 'site']).set_index(['subj', 'strand'], drop=True)

    adapter = {}
    logger.info(" - reading adapter fasta ...")
    for record in SeqIO.parse(args.adapter, "fasta"):
        adapter[str(record.id)] = str(record.seq)
    adapter = pd.DataFrame.from_dict(adapter, orient='index', columns=['seq'], dtype='string')
    logger.info('{:>11} adapter sequence(s) in fasta file'.format(len(adapter)))

    logger.info(" - reading genome fasta ...")
    genome = {}
    for record in SeqIO.parse(args.genome, "fasta"):
        genome[str(record.id)] = str(record.seq)
    genome = pd.DataFrame.from_dict(genome, orient='index', columns=['seq'], dtype='string')
    logger.info('{:>11} genomic sequence(s) in fasta file'.format(len(genome)))

    logger.info(" - reading fastq files ...")
    reads = {}
    for fqFile in fq_files:
        for record in SeqIO.parse(fqFile, "fastq"):
            reads[str(record.id)] = str(record.seq)
    reads = pd.DataFrame.from_dict(reads, orient='index', columns=['seq'], dtype='string')
    logger.info("{:>11} reads in dataset\n".format(len(reads)))

    if not args.adapter_alignment:
        logger.info(" - performing reads to adapter reference mapping ...")
        fq_fn = " ".join(fq_files)
        ref_fn = args.adapter
        args.adapter_alignment = "{}.adapter_alignment.paf".format(args.prefix)
        exit_code = run_minimap2(ref_fn, fq_fn, args.adapter_alignment)
        if exit_code:
            logger.error('ERROR: adapter reference mapping failed with exit code', exit_code)
            exit(1)
    if args.adapter_alignment.endswith('.paf'):
        ad_algn_df = parse_paf(args.adapter_alignment, cigar=True)#.set_index('qid')
    elif args.adapter_alignment.endswith('.bam'):
        ad_algn_df = parse_bam(args.adapter_alignment, reads)#.set_index('qid')
    else:
        logger.error('Alignment file type not supported: {}'.format(args.adapter_alignment))
        exit(1)
    breakpoint()
    logger.info("{:>11} {:>7} primary alignments against adapter sequence(s)".format(len(ad_algn_df), ""))
    c = sum(ad_algn_df.strand == '+')
    logger.info("{:>11} {:>5.1f} % against (+) strand".format(c, c/len(ad_algn_df)*100.))
    c = sum(ad_algn_df.strand == '-')
    logger.info("{:>11} {:>5.1f} % against (-) strand".format(c, c/len(ad_algn_df)*100.))
    c = len(set(ad_algn_df.qid))
    logger.info("{:>11} {:>5.1f} % of reads align against any adapter sequence".format(c, c/len(reads)*100.))

    if not args.genome_alignment:
        logger.info(" - performing reads to genome reference mapping ...")
        ref_fn = args.genome
        args.genome_alignment = "{}.genome_alignment.paf".format(args.prefix)
        exit_code = run_minimap2(ref_fn, fq_fn, args.genome_alignment)
        if exit_code:
            logger.error('ERROR: adapter reference mapping failed with exit code', exit_code)
            exit(1)
    if args.genome_alignment.endswith('.paf'):
        gn_algn_df = parse_paf(args.genome_alignment, cigar=True)#.set_index('qid')
    elif args.genome_alignment.endswith('.bam'):
        gn_algn_df = parse_bam(args.genome_alignment, reads) #.set_index('qid')
    else:
        logger.error('Alignment file type not supported: {}'.format(args.genome_alignment))
        exit(1)
    logger.info("{:>11} {:>7} primary alignments against genomic sequence(s)".format(len(gn_algn_df), ""))
    c = sum(gn_algn_df.strand == '+')
    logger.info("{:>11} {:>5.1f} % against (+) strand".format(c, c/len(gn_algn_df)*100.))
    c = sum(gn_algn_df.strand == '-')
    logger.info("{:>11} {:>5.1f} % against (-) strand".format(c, c/len(gn_algn_df)*100.))
    c = len(set(gn_algn_df.qid))
    logger.info("{:>11} {:>5.1f} % of reads align against any genome sequence".format(c, c/len(reads)*100.))

    if not (args.mean and args.std):
        logger.info(" - determining statistics about per-base difference ([align. subject len] - [align. query len]) / [align. query len] from genome alignments")
        sequence_length_stats(gn_algn_df)
        logger.info("{:>11.4f} mean of per-base difference in sequence length".format(args.mean))
        logger.info("{:>11.4f} std dev of per-base difference in sequence length".format(args.std))

    logger.info(" - joining alignment data and filtering for adjacent adapter-genome alignments")
    # determine order of alignments of each read with respect to each adapter-subject pair
    d_ad = pd.merge(pd.DataFrame(reads.index, columns=['rid']).set_index('rid', drop=False), 
                    ad_algn_df.set_index('qid'),
                    how='outer', left_index=True, right_index=True, sort=False)
    df = pd.merge(d_ad, 
                  gn_algn_df.set_index('qid'), 
                  how='outer', left_index=True, right_index=True, sort=False, suffixes=('_ad', '_gn'))
    # reset index
    df = df.reset_index(drop=True)
    sel = df.subj_ad.notnull() & \
          df.subj_gn.notnull()

    logger.info('{:>11} {:>7} entities after joining'.format(len(df), ""))
    c = sum(np.logical_not(sel))
    logger.info('{:>11} {:>5.1f} % not aligning against both an adapter and a genomic sequence'.format(c, c/len(df)*100.))
    c = len(set(df.loc[sel, 'rid']))
    logger.info('{:>11} {:>5.1f} % reads remaining that contain potential transitions between adapter and genomic sequence.'.format(c, c/len(reads)*100.))

    # identify rows describing a transition from adapter to genomic sequence or vise versa
    df.loc[sel,'trans_order'] = 0 # qst_ad == qst_gn
    df.loc[sel & (df.qst_ad < df.qst_gn), 'trans_order'] = 1 # adapter -> genome
    df.loc[sel & (df.qst_ad > df.qst_gn), 'trans_order'] = -1 # genome -> adapter

    # remove pairs of alignments which are too distant from one another
    sel_ = (sel & (df.trans_order ==  1.) & ((df.qst_gn - df.qen_ad) <= args.max_dist)) | \
           (sel & (df.trans_order == -1.) & ((df.qst_ad - df.qen_gn) <= args.max_dist))
    c = sum(sel & np.logical_not(sel_))
    logger.info('{:>11} {:>5.1f} % transitions have local alignments that are >{} nt apart from each other'.format(c, c/len(df)*100., args.max_dist))
    sel = sel_
    c = sum(sel)
    logger.info('{:>11} {:>5.1f} % potential transitions remaining'.format(c, c/len(df)*100.))

    # deleting non-potential transitions at this point so that all procentual values are with respect to potential ones only
    logger.info(" - deleting all entries from dataframe that are not potential transitions")
    dtypes = {"{}_{}".format(key, suffix):np.int32 for suffix in ['ad', 'gn'] \
              for key in ["qlen","qst","qen","slen","sst","sen","mlen","blen","mapq"]}
    df = df.drop(df.index[np.logical_not(sel)]).astype(dtypes)
    sel = df.subj_ad.notnull() & \
          df.subj_gn.notnull()

    
    c = sum(df[sel].qst_ad < df[sel].qst_gn)
    logger.info('{:>11} {:>5.1f} % of total potential transpositions are adapter -> genome'.format(c, c/len(df[sel])*100.))
    c = sum(df[sel].qst_ad > df[sel].qst_gn)
    logger.info('{:>11} {:>5.1f} % of total potential transpositions are genome -> adapter'.format(c, c/len(df[sel])*100.))

    # select all entries with a sufficiently long alignment to both the adapter and the genome
    logger.info(' - filter potential subject transitions based on alignment lengths')
    sel_ = (df.blen_ad >= args.min_adapter_blen) & \
           (df.blen_gn >= args.min_genome_blen)
    c = sum(sel & np.logical_not(sel_))
    logger.info('{:>11} {:>5.1f} % of total potential transpositions filtered'.format(c, c/len(df)*100.))
    c = sum(sel & np.logical_not(df.blen_ad >= args.min_adapter_blen))
    logger.info('{:>11} {:>5.1f} % bc adapter alignment length < {}'.format(c, c/len(df)*100., args.min_adapter_blen))
    c = sum(sel & np.logical_not(df.blen_gn >= args.min_genome_blen))
    logger.info('{:>11} {:>5.1f} % bc genome alignment length < {}'.format(c, c/len(df)*100., args.min_genome_blen))
    sel &= sel_
    c = sum(sel)
    logger.info('{:>11} {:>5.1f} % potential transitions remaining'.format(c, c/len(df)*100.))

    # determine query seq start and end in read coordinates
    logger.info(" - determine query seq start and end in read coordinates")
    df = set_qst_and_qen(df, sel)
    c = sum(sel & (df.qst.isnull() | df.qen.isnull()))
    logger.info('{:>11} {:>5.1f} % of remaining potential transpositions had < {} matching terminal bases'.format(c, c/sum(sel)*100., args.wordsize))
    sel &= df.qst.notnull() & df.qen.notnull()
    c = sum(sel)
    logger.info('{:>11} {:>5.1f} % potential transitions remaining'.format(c, c/len(df)*100.))

    # filter out entries with query seq. that are too long -> distance between adapter and genome alignment is too large
    logger.info(" - exclude entries based on query sequence length (adapter-genome alignment distance) from analysis")
    df.loc[sel, 'qlen'] = df[sel].qen - df[sel].qst
    too_short = df.qlen < 0.
    too_long = df.qlen > args.max_dist
    c = sum(too_short)
    logger.info('{:>11} {:>5.1f} % of remaining potential transpositions removed due to query seq. length < 0 nt (alignment overlap of more than {} nt)'.format(c, c/sum(sel)*100., args.strip))
    c = sum(too_long)
    logger.info('{:>11} {:>5.1f} % of remaining potential transpositions removed due to query seq. length > {} nt'.format(c, c/sum(sel)*100., args.max_dist))
    sel &= np.logical_not(too_short) & np.logical_not(too_long)
    c = sum(sel)
    logger.info('{:>11} {:>5.1f} % potential transitions remaining'.format(c, c/len(df)*100.))

    ## set query seq for each row
    #logger.info(" - determining query sequences")
    #df.loc[sel, 'query_seq'] = df[sel].parallel_apply(lambda row: reads.loc[row.rid].seq[int(row.qst) : int(row.qen)], axis=1)
    #
    logger.info(" - determining reference sequences")
    df.loc[sel, 'max_ref_len'] = ((df[sel].qen - df[sel].qst) * (1 + args.mean + 3 * args.std) + 0.5).round()
    df = set_references(df, sel)
    #df.loc[sel, 'ref1'] = df.loc[sel].parallel_apply(lambda row: get_ref1(row), axis=1)
    #df.loc[sel, 'ref2'] = df.loc[sel].parallel_apply(lambda row: get_ref2(row), axis=1)

    #if args.sites_of_interest:
    #    logger.info(" - reject potential transitions that do not span a site of interest")
    #    soi = pd.read_csv(args.sites_of_interest, header=None, sep='\t', names=['subj', 'site'])
    #    #sel_ = pd.Series(np.full(shape=sel.shape, fill_value=False))
    #    #for i,(subj, site) in tqdm(soi.iterrows()):
    #    #    sel_ref1 = (df.subj_ref1 == subj) & (df.sst_ref1 <= site) & (df.sen_ref1 > site)
    #    #    sel_ref2 = (df.subj_ref2 == subj) & (df.sst_ref2 <= site) & (df.sen_ref2 > site)
    #    #    sel_ |= sel_ref1 | sel_ref2
    #    c = sum(sel_)
    #    logger.info('{:>11} {:>5.1f} % of remaining potential transpositions rejected'.format(c, c/sum(sel)*100.))
    #    sel &= sel_
    #    c = sum(sel)
    #    logger.info('{:>11} {:>5.1f} % potential transitions remaining'.format(c, c/len(df)*100.))

    logger.info(' - aligning')
    ct.CDLL(clib).init(args.sequence_type == 'nucl')

    # create arrays that are passed to c functions
    max_subj_len = 2*int((args.max_dist * (1 + args.mean + 3 * args.std) + 0.5).round())

    scores = np.empty(shape=(max_subj_len+1, args.max_dist+1), dtype=np.int16)
    ops = np.empty(shape=(6,max_subj_len+1, args.max_dist+1), dtype=np.bool_)
    transitions = np.empty(shape=(2, max_subj_len+1, max_subj_len+1), dtype=np.int16)
    align_ends = np.empty(shape=(2, max_subj_len+1, args.max_dist+1), dtype=np.bool_)
    reachable = np.empty(shape=(2*max_subj_len+1, args.max_dist+1), dtype=np.bool_)
    cigarbuffer = ct.create_string_buffer(args.max_dist + 2 * max_subj_len + 1)
    for arr in [scores, ops, transitions, align_ends, reachable]:
        if arr.flags['C_CONTIGUOUS'] == False:
            arr = np.ascontiguousarray(arr, dtype=arr.dtype)

    scores_pp = (scores.ctypes.data + np.arange(scores.shape[0]) * scores.strides[0]).astype(np.uintp)
    ops_pp = [(ops[i].ctypes.data + np.arange(ops[i].shape[0]) * ops[i].strides[0]).astype(np.uintp) for i in range(6)]
    transitions_pp0 = (transitions[0].ctypes.data + np.arange(transitions[0].shape[0]) * transitions[0].strides[0]).astype(np.uintp)
    transitions_pp1 = (transitions[1].ctypes.data + np.arange(transitions[1].shape[0]) * transitions[1].strides[0]).astype(np.uintp)
    align_ends_pp0 = (align_ends[0].ctypes.data + np.arange(align_ends[0].shape[0]) * align_ends[0].strides[0]).astype(np.uintp)
    align_ends_pp1 = (align_ends[1].ctypes.data + np.arange(align_ends[1].shape[0]) * align_ends[1].strides[0]).astype(np.uintp)
    reachable_pp = (reachable.ctypes.data + np.arange(reachable.shape[0]) * reachable.strides[0]).astype(np.uintp)
    
    # initialize the fields that never change
    scores[0,0] = 0
    ops[:,0,:] = False
    ops[:,:,0] = False
    ops[0,0,1:] = True
    ops[2,0,0] = True
    ops[5,0,0] = True
    if args.gap_extension == 0:
        scores[0,1:] = np.full(args.max_dist, args.gap_open)
    else:
        scores[0,1:] = np.arange(args.gap_open + args.gap_extension, 
                                 args.gap_open + (args.max_dist+1)*args.gap_extension, 
                                 args.gap_extension)

    # start the alignment
    df = df.parallel_apply(lambda row: c_align_row(row), axis=1)

    sel = sel & df.score.notnull()
    c = sum(sel)
    logger.info('{:>11} {:>5.1f} % of potential transitions successfully aligned'.format(c, c/len(df)*100.))

    df.loc[sel & (df.qlen > 0), 'norm_score'] = df[sel].score / df[sel].qlen
    df.loc[sel & (df.qlen == 0), 'norm_score'] = args.match # if the two alignments are fitting together perfetly

    for upper, lower in [(1., .9), (.9, .8), (.8, .7), (.7, .6), (.6, .5), (.5, 0.)]:
        c = sum(sel & (upper >= df.norm_score) & (df.norm_score > lower))
        logger.info('{:>11} {:>5.1f} % with {} >= normed score > {}'.format(c, c/sum(sel)*100., upper, lower))

    if args.plot:
        plot_norm_score_distribution(df[sel], "all data")

    fn = args.prefix + ".alignment.df.pkl"
    logger.info(' - writing results to file {}'.format(fn))
    df.loc[sel].to_pickle(fn)

