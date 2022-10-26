from crossalign import *
import argparse, os, re
from matplotlib import pyplot as plt
import logging
import Bio

def parse_args(args_=None):
    global args
    parser = argparse.ArgumentParser(description='Outputs per feature hits based on a genbank annotation file and a crossalign output.',
                                     formatter_class=ArgHelpFormatter, 
                                     add_help=False)

    main_group = parser.add_argument_group('Main Arguments')
    main_group.add_argument('--results',
                            required=True,
                            help='''Path to a .pkl file output by crossalign. If multiple are given, 
                            their output is combined. In this case, a new prefix must be given.''')
    main_group.add_argument('--genbank',
                            required=True,
                            help='''Path to a genbank file containing the genome annotation.''')
    main_group.add_argument('--min_norm_score',
                            default=np.NINF,
                            help='''Reject transitions with a length normalized alignment score lower than this threshold.''')
    main_group.add_argument('--filter_ambiguous',
                            action="store_true",
                            help='''Reject transitions that are ambiguous with respect to the genomic insertion site.''')
    main_group.add_argument('--promotor_bases',
                            type=int,
                            default=50,
                            help='''Number of included bases upstream of the feature start.''')
    main_group.add_argument('--essential_fraction',
                            type=float,
                            default=0.9,
                            help='''Fraction of feature length to be considered essential.''')
    main_group.add_argument('--rejected_genbank_keys',
                            nargs='*',
                            default=['source', 'gene', 'rRNA'],
                            help='''Reject Genbank features with one of the listed key types.''')

    plot_group = parser.add_argument_group('Plotting Arguments')
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
    plot_group.add_argument('--bins',
                            default=60,
                            help='number of bins in the distribution plot.')

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
    args.prefix = re.fullmatch('(.*)\.alignment\.df\.pkl', os.path.basename(args.results)).group(1)

    return args

def gb_to_feature_list(fp, exclude_keys=['source', 'gene', 'rRNA'], prom_dist=50, essential=0.8):
    data = []
    with open(fp) as handle:
        for record in Bio.GenBank.parse(handle):
            replicon = record.version
            for feature in record.features:
                if feature.key == "source":
                    repl_len = int(re.match("(\d+)\.\.(\d+)", feature.location).group(2))
                if feature.key not in exclude_keys:
                    locus_tag = [qf.value for qf in feature.qualifiers if qf.key == '/locus_tag=']
                    if not locus_tag:
                        locus_tag = "_".join([feature.key, feature.location])
                        if locus_tag in [locus_tag for replicon, locus_tag, gene, start, end, strand, repl_len in data]:
                            i = 2
                            while f"{locus_tag}_{i}" in [locus_tag for replicon, locus_tag, gene, start, end, strand, repl_len in data]:
                                i += 1
                            locus_tag = f"{locus_tag}_{i}"
                        print("WARNING: no locus_tag for feature", feature.key, feature.location, " --> set to", locus_tag)
                    else:
                        locus_tag = locus_tag[0].strip('"')
                    gene = [qf.value for qf in feature.qualifiers if qf.key == '/gene=']
                    gene = gene[0].strip('"') if gene else ""
                    location = re.match("(complement\()?(join\()?\<?(\d+)\.\.(\d+,\d+\.\.)*\>?(\d+)\)*", feature.location)
                    start = int(location.group(3))
                    end = int(location.group(5))
                    strand = "-" if feature.location.startswith('complement') else '+'
                    data.append((replicon, locus_tag, gene, start, end, strand, repl_len))
    df_ = pd.DataFrame(data, columns=['replicon', 'locus_tag', 'gene', 'start', 'end', 'strand', 'repl_len']).set_index(['replicon', 'locus_tag'])
    # temporarily change start and stop of features that span the ORI
    df_.loc[df_.end < df_.start, 'end'] = df_.loc[df_.end < df_.start, 'repl_len'] + df_.loc[df_.end < df_.start, 'end']
    # set promotor region start
    df_.loc[df_.strand == '+', 'prom'] = df_.loc[df_.strand == '+', 'start'] - prom_dist
    df_.loc[df_.strand == '-', 'prom'] = df_.loc[df_.strand == '-', 'end'] + prom_dist
    df_.loc[df_.prom < 1, 'prom'] = df_.loc[df_.prom < 1, 'repl_len'] + df_.loc[df_.prom < 1, 'prom']
    df_.loc[df_.prom > df_.repl_len, 'prom'] = - df_.loc[df_.prom > df_.repl_len, 'prom']
    # set essential region end
    df_.loc[df_.strand == '+', 'essen'] = df_.loc[df_.strand == '+', 'start'] +\
        (essential * (df_.loc[df_.strand == '+', 'end'] - df_.loc[df_.strand == '+', 'start'])).round()
    df_.loc[df_.strand == '-', 'essen'] = df_.loc[df_.strand == '-', 'end'] -\
        (essential * (df_.loc[df_.strand == '-', 'end'] - df_.loc[df_.strand == '-', 'start'])).round()
    df_.loc[df_.essen < 1, 'essen'] = df_.loc[df_.essen < 1, 'repl_len'] + df_.loc[df_.essen < 1, 'essen']
    df_.loc[df_.essen > df_.repl_len, 'essen'] = - df_.loc[df_.essen > df_.repl_len, 'essen']
    # redo changes to start and stop of features that span the ORI
    df_.loc[df_.end > df_.repl_len, 'end'] = df_.loc[df_.end > df_.repl_len, 'end'] - df_.loc[df_.end > df_.repl_len, 'repl_len']
    # set new start and end positions for the hit region
    df_.loc[df_.strand == '+', "reg_start"] = df_.loc[df_.strand == '+', "prom"]
    df_.loc[df_.strand == '+', "reg_end"] = df_.loc[df_.strand == '+', "essen"]
    df_.loc[df_.strand == '-', "reg_start"] = df_.loc[df_.strand == '-', "essen"]
    df_.loc[df_.strand == '-', "reg_end"] = df_.loc[df_.strand == '-', "prom"]
    # compute region length
    df_.loc[df_.reg_end >= df_.reg_start, 'reg_len'] = df_.loc[df_.reg_end >= df_.reg_start, 'reg_end'] - df_.loc[df_.reg_end >= df_.reg_start, 'reg_start'] + 1
    df_.loc[df_.reg_end < df_.reg_start, 'reg_len'] = df_.loc[df_.reg_end < df_.reg_start, 'reg_end'] + (df_.loc[df_.reg_end < df_.reg_start, 'repl_len'] - df_.loc[df_.reg_end < df_.reg_start, 'reg_start']) + 1
    return df_.drop(['prom', 'essen'], axis=1)

def count_feature_hits(features, d, min_norm_score=np.NINF, filter_ambiguous=False):
    hits = {}
    sum_hits = {}
    for subj in d.subj_gn.unique():
        hits[subj] = {'+' : np.zeros(features.loc[subj].iloc[0].repl_len + 1),
                      '-' : np.zeros(features.loc[subj].iloc[0].repl_len + 1)}
        sum_hits[subj] = 0
        sel = (d.subj_gn == subj) & (d.norm_score >= min_norm_score)
        if filter_ambiguous:
            sel &= (d.amb == False)
        for strand in '+-':
            sel_ = sel & (d.strand_gn == strand)
            d_ = d.loc[sel_].groupby('site_gn')[['one_div_sites']].sum()
            hits[subj][strand][d_.index.astype(int)] = d_.one_div_sites
            sum_hits[subj] += hits[subj][strand].sum()
    
    def region_hits(row):
        start, end = int(row.reg_start), int(row.reg_end)
        if start <= end:
            row['hits_+'] = hits[row.name[0]]['+'][start-1:end].sum()
            row['hits_-'] = hits[row.name[0]]['-'][start-1:end].sum()
        else:
            row['hits_+'] = hits[row.name[0]]['+'][start-1:].sum() + hits[row.name[0]]['+'][:end].sum()
            row['hits_-'] = hits[row.name[0]]['-'][start-1:].sum() + hits[row.name[0]]['-'][:end].sum()
        row['hits'] = row['hits_+'] + row['hits_-']
        row['hits_norm'] = row.hits / (sum_hits[row.name[0]] * row.reg_len)
        return row
    
    return features.apply(region_hits, axis=1)

def main():
    # load crossalign results
    df = pd.read_pickle(args.results)
    df['one_div_sites'] = 1 / df.transitions.str.len()
    d = df.explode('transitions')

    d['site_gn'] = np.nan
    d['site_ad'] = np.nan
    d.loc[(d.order == 1), 'site_gn'] = d.loc[(d.order == 1)].transitions.str[1]
    d.loc[(d.order == -1), 'site_gn'] = d.loc[(d.order == -1)].transitions.str[0]
    d.loc[(d.order == 1), 'site_ad'] = d.loc[(d.order == 1)].transitions.str[0]
    d.loc[(d.order == -1), 'site_ad'] = d.loc[(d.order == -1)].transitions.str[1]

    features = gb_to_feature_list(args.genbank, prom_dist=args.promotor_bases, essential=args.essential_fraction)
    hits = count_feature_hits(features, d, min_norm_score=args.min_norm_score, filter_ambiguous=args.filter_ambiguous)

    fp = f"{args.prefix}.feat_hits.csv"
    print(f"writing to file {fp} ...")
    hits.to_csv(fp)

    fig,ax = plt.subplots(figsize=(args.width, args.height))
    ax.hist(hits.hits_norm, bins=args.bins)
    fp = f"{args.prefix}.norm_feat_hits_dist.{args.file_type}"
    print(f"writing to file {fp} ...")
    plt.savefig(fp, dpi=args.dpi)

if __name__ == "__main__":
    args = parse_args()
    main()
