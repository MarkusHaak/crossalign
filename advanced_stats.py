from crossalign import *
import argparse, os, re
from matplotlib import pyplot as plt
import logging

def parse_args(args_=None):
    global args
    parser = argparse.ArgumentParser(description='Creates additional stats and plots for crossalign output files.',
                                     formatter_class=ArgHelpFormatter, 
                                     add_help=False)

    main_group = parser.add_argument_group('Main Arguments')
    main_group.add_argument('--reads',
                            required=True,
                            nargs='+',
                            help='fastq file(s) or path to directory(s) containing fastq files (recursion depth 1).')
    main_group.add_argument('--genome',
                            required=True,
                            help='First set of reference sequences (name deprecated, can be any (multiple) fasta file.')
    main_group.add_argument('--adapter',
                            required=True,
                            help='Second set of reference sequences (name deprecated, can be any (multiple) fasta file.')
    main_group.add_argument('--results',
                            required=True,
                            help='path to the .pkl file output by crossalign.')
    
    plot_group = parser.add_argument_group('Plotting Arguments')
    plot_group.add_argument('--outdir',
                            help='directory where the output plots are saved. (default: input directory)')
    plot_group.add_argument('--width',
                            type=float,
                            default=6.0,
                            help='Width of output plots in inches.')
    plot_group.add_argument('--height',
                            type=float,
                            default=4.0,
                            help='Height of output plots in inches.')
    plot_group.add_argument('--file_type',
                            choices=['png', 'jpg', 'svg', 'pdf'],
                            default='png',
                            help='Output format of generated plots.')
    plot_group.add_argument('--dpi',
                            type=int,
                            default=150,
                            help='Resolution in dpi.')
    plot_group.add_argument('--lower_quantile',
                            default=0.01,
                            help='Only data between the lower and upper quantile will be plotted.')
    plot_group.add_argument('--upper_quantile',
                            default=0.99,
                            help='Only data between the lower and upper quantile will be plotted.')
    plot_group.add_argument('--read_length_bins',
                            default=60,
                            help='number of bins in the read-length distribution plot.')

    help_group = parser.add_argument_group('Help')
    help_group.add_argument('-h', '--help', 
                            action='help', 
                            default=argparse.SUPPRESS,
                            help='Show this help message and exit.')
    if args_:
        args = parser.parse_args(args_)
    else:
        args = parser.parse_args()

    if not args.results.endswith('.alignment.df.pkl') or not os.path.exists(args.results):
        print("Please select a .alignment.df.pkl file as the --results file.")
        exit(1)
    if not args.outdir:
        args.outdir = os.path.dirname(os.path.abspath(args.results))
    else:
        args.outdir = os.path.abspath(args.outdir)
    args.prefix = re.fullmatch('(.*)\.alignment\.df\.pkl', os.path.basename(args.results)).group(1)

    if os.path.exists(args.outdir):
        if not os.path.isdir(args.outdir):
            print(f"Path does exist but is not a directory: {args.outdir}")
            exit(1)
        elif not os.access(args.outdir, os.W_OK):
            print(f"Directory is not writable: {args.outdir}")
            exit(1)
    else:
        os.makedirs(args.outdir)
    return args

def quantil_barplot(d, x_lbl, filename, low=0.01, high=0.99, bins=50):
    '''Creates a barplot from the data between the lower and upper quantil in a pandas Series d'''
    low_, high_ = d.quantile([low,high])
    data = d.loc[(d >= low_) & (d <= high_)]
    fig, ax = plt.subplots(figsize=(args.width, args.height))
    ax.grid(alpha=0.2)
    ax.hist(data, bins=bins)
    ax.set(xlabel=x_lbl)
    fp = os.path.join(args.outdir, f"{args.prefix}.{filename}.{args.file_type}")
    plt.savefig(fp, dpi=args.dpi)


def main():
    # load sequence data
    adapter, genome, reads = read_sequence_data(args.reads, args.adapter, args.genome)
    # load crossalign results
    df = pd.read_pickle(args.results)
    d = df.explode('transitions')
    # add site columns
    d['site_gn'] = np.nan
    d['site_ad'] = np.nan
    d.loc[(d.order == 1), 'site_gn'] = d.loc[(d.order == 1)].transitions.str[1]
    d.loc[(d.order == -1), 'site_gn'] = d.loc[(d.order == -1)].transitions.str[0]
    d.loc[(d.order == 1), 'site_ad'] = d.loc[(d.order == 1)].transitions.str[0]
    d.loc[(d.order == -1), 'site_ad'] = d.loc[(d.order == -1)].transitions.str[1]
    # print stats to std:
    print(f"{sum([len(reads[i]) for i in reads]):>10} total reads")
    print(f"{df.rid.unique().size:>10} mapped reads")
    print(f"{len(df.loc[df.transitions.str.len() > 1]):>10} mapped reads with multiple transitions")
    print(f"{df.index.unique().size:>10} total transitions")
    print(f"{len(df.loc[df.amb == False]):>10} unambiguous transitions (with respect to the transition site)")
    print(f"{len(d):>10} transition site combinations (all transitions)")
    print(f"{len(d.loc[d.amb == False]):>10} transition site combinations (unique transitions)")
    print()
    print(f'{d.set_index(["strand_ad","site_ad","subj_ad"]).index.unique().size:>10} unique sites on any adapter sequence (all transitions)')
    for subj_ad in adapter.index:
        d_sel = d.loc[(d.subj_ad == subj_ad)]
        print(f'{d_sel.set_index(["strand_ad","site_ad"]).index.unique().size:>10} unique sites on adapter sequence "{subj_ad}" (all transitions):')
        print(f'{d_sel.loc[d_sel.strand_ad == "+"].site_ad.unique().size:>10}     on (+) strand')
        print(f'{d_sel.loc[d_sel.strand_ad == "-"].site_ad.unique().size:>10}     on (-) strand')
    print(f'{d.loc[d.amb == False].set_index(["strand_ad","site_ad","subj_ad"]).index.unique().size:>10} unique sites on any adapter sequence (unambiguous transitions)')
    for subj_ad in adapter.index:
        print(f'{d_sel.loc[(d_sel.amb == False)].set_index(["strand_ad","site_ad"]).index.unique().size:>10} unique sites on adapter sequence "{subj_ad}" (unambiguous transitions)')
        print(f'{d_sel.loc[(d_sel.amb == False) & (d_sel.strand_ad == "+")].site_ad.unique().size:>10}     on (+) strand')
        print(f'{d_sel.loc[(d_sel.amb == False) & (d_sel.strand_ad == "-")].site_ad.unique().size:>10}     on (-) strand')
    print()
    print(f'{d.set_index(["strand_gn","site_gn","subj_gn"]).index.unique().size:>10} unique sites on any genome sequence (all transitions)')
    for subj_gn in genome.index:
        d_sel = d.loc[(d.subj_gn == subj_gn)]
        print(f'{d_sel.set_index(["strand_gn","site_gn"]).index.unique().size:>10} unique sites on genome sequence "{subj_gn}" (all transitions):')
        print(f'{d_sel.loc[d_sel.strand_gn == "+"].site_gn.unique().size:>10}     on (+) strand')
        print(f'{d_sel.loc[d_sel.strand_gn == "-"].site_gn.unique().size:>10}     on (-) strand')
    print(f'{d.loc[d.amb == False].set_index(["strand_gn","site_gn","subj_gn"]).index.unique().size:>10} unique sites on any genome sequence (unambiguous transitions)')
    for subj_gn in genome.index:
        d_sel = d.loc[(d.subj_gn == subj_gn)]
        print(f'{d_sel.loc[(d_sel.amb == False)].set_index(["strand_gn","site_gn"]).index.unique().size:>10} unique sites on genome sequence "{subj_gn}" (unambiguous transitions)')
        print(f'{d_sel.loc[(d_sel.amb == False) & (d_sel.strand_gn == "+")].site_gn.unique().size:>10}     on (+) strand')
        print(f'{d_sel.loc[(d_sel.amb == False) & (d_sel.strand_gn == "-")].site_gn.unique().size:>10}     on (-) strand')
    
    # create and save some plots:
    quantil_barplot(reads.seq.str.len(), "read length / nt", "read_length_dist", bins=args.read_length_bins)
    quantil_barplot(df.norm_score, "length normalized score", "length_norm_score_dist", high=1.)
    quantil_barplot(df.transitions.str.len(), "transitions per read", "transitions_per_read_dist", 
                    low=0., high=1., bins=max(1, df.transitions.str.len().max()-1))
    quantil_barplot(d.groupby(["strand_gn","site_gn"]).size(), "transitions per site (all)", "transitions_per_site_dist_all", 
                    low=0., bins=max(1, int(d.groupby(["strand_gn","site_gn"]).size().quantile(args.lower_quantile)-1)))
    quantil_barplot(d.loc[(d.amb == False) & (d.norm_score >= 0.0)].groupby(["strand_gn","site_gn"]).size(), "transitions per site (all)", "transitions_per_site_dist_uniq", 
                    low=0., bins=max(1, int(d.loc[(d.amb == False) & (d.norm_score >= 0.0)].groupby(["strand_gn","site_gn"]).size().quantile(args.lower_quantile)-1)))

if __name__ == "__main__":
    args = parse_args()
    main()