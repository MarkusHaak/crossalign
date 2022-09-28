from crossalign import *
import argparse, os, re
from matplotlib import pyplot as plt
import logging

def parse_args(args_=None):
    global args
    parser = argparse.ArgumentParser(description='Creates wiggle files from crossalign output files.',
                                     formatter_class=ArgHelpFormatter, 
                                     add_help=False)

    main_group = parser.add_argument_group('Main Arguments')
    main_group.add_argument('--results',
                            nargs='+',
                            help='''Path(s) to .pkl file(s) output by crossalign. If multiple are given, 
                            their output is combined. In this case, a new prefix must be given.''')
    main_group.add_argument('--prefix',
                            help='Prefix for the created output files. (Default: same as prefix of single input file)')

    help_group = parser.add_argument_group('Help')
    help_group.add_argument('-h', '--help', 
                            action='help', 
                            default=argparse.SUPPRESS,
                            help='Show this help message and exit.')
    if args_:
        args = parser.parse_args(args_)
    else:
        args = parser.parse_args()

    for fp in args.results:
        if not fp.endswith('.alignment.df.pkl') or not os.path.exists(fp):
            print("Please select a .alignment.df.pkl file as the --results file.")
            exit(1)
    if len(args.results) == 1 and not args.prefix:
        args.prefix = re.fullmatch('(.*)\.alignment\.df\.pkl', args.results[0]).group(1)
    elif len(args.results) > 1 and not args.prefix:
        print("Please specify a new --prefix if using multiple results input files.")
        exit(1)

    return args

def main():
    dfs, ds = [], []
    for results in args.results:
        # load crossalign results
        df = pd.read_pickle(results)
        df['one_div_sites'] = 1 / df.transitions.str.len()
        d = df.explode('transitions')
        # add site columns
        d['site_gn'] = np.nan
        d['site_ad'] = np.nan
        d.loc[(d.order == 1), 'site_gn'] = d.loc[(d.order == 1)].transitions.str[1]
        d.loc[(d.order == -1), 'site_gn'] = d.loc[(d.order == -1)].transitions.str[0]
        d.loc[(d.order == 1), 'site_ad'] = d.loc[(d.order == 1)].transitions.str[0]
        d.loc[(d.order == -1), 'site_ad'] = d.loc[(d.order == -1)].transitions.str[1]
        dfs.append(df)
        ds.append(d)
    
    # all
    # genome by subj and strand
    sites = ds[0].groupby(["subj_gn", "strand_gn", "site_gn"])[['one_div_sites']].sum()
    for d in ds[1:]:
        sites = sites.add(d.groupby(["subj_gn", "strand_gn", "site_gn"])[['one_div_sites']].sum(),
                          fill_value=0.)
    for subj_gn in d.subj_gn.unique():
        for strand in ['+', '-']:
            fp = f"{args.prefix}.all.genome.{subj_gn}_{strand}.wig"
            print(f"writing to file {fp} ...")
            with open(fp, 'w') as f:
                print(f"variableStep  chrom={subj_gn}", file=f)
                for site, row in sites.loc[subj_gn, strand].iterrows():
                    print(int(site)+1, row.one_div_sites, file=f)
    # genome by subj and strands combined
    sites = ds[0].groupby(["subj_gn", "site_gn"])[['one_div_sites']].sum()
    for d in ds[1:]:
        sites = sites.add(d.groupby(["subj_gn", "site_gn"])[['one_div_sites']].sum(),
                          fill_value=0.)
    for subj_gn in d.subj_gn.unique():
        fp = f"{args.prefix}.all.genome.{subj_gn}.wig"
        print(f"writing to file {fp} ...")
        with open(fp, 'w') as f:
            print(f"variableStep  chrom={subj_gn}", file=f)
            for site, row in sites.loc[subj_gn].iterrows():
                print(int(site)+1, row.one_div_sites, file=f)
    # adapter by subj and strand
    sites = ds[0].groupby(["subj_ad", "strand_ad", "site_ad"])[['one_div_sites']].sum()
    for d in ds[1:]:
        sites = sites.add(d.groupby(["subj_ad", "strand_ad", "site_ad"])[['one_div_sites']].sum(),
                          fill_value=0.)
    for subj_ad in d.subj_ad.unique():
        for strand in ['+', '-']:
            fp = f"{args.prefix}.all.adapter.{subj_ad}_{strand}.wig"
            print(f"writing to file {fp} ...")
            with open(fp, 'w') as f:
                print(f"variableStep  chrom={subj_ad}", file=f)
                for site, row in sites.loc[subj_ad, strand].iterrows():
                    print(int(site)+1, row.one_div_sites, file=f)
    # adapter by subj and strands combined
    sites = ds[0].groupby(["subj_ad", "site_ad"])[['one_div_sites']].sum()
    for d in ds[1:]:
        sites = sites.add(d.groupby(["subj_ad", "site_ad"])[['one_div_sites']].sum(),
                          fill_value=0.)
    for subj_ad in d.subj_ad.unique():
        fp = f"{args.prefix}.all.adapter.{subj_ad}.wig"
        print(f"writing to file {fp} ...")
        with open(fp, 'w') as f:
            print(f"variableStep  chrom={subj_ad}", file=f)
            for site, row in sites.loc[subj_ad].iterrows():
                print(int(site)+1, row.one_div_sites, file=f)

    # unambiguous
    # genome by subj and strand
    sites = ds[0].loc[d.amb == False].groupby(["subj_gn", "strand_gn", "site_gn"])[['one_div_sites']].sum()
    for d in ds[1:]:
        sites = sites.add(d.loc[d.amb == False].groupby(["subj_gn", "strand_gn", "site_gn"])[['one_div_sites']].sum(),
                          fill_value=0.)
    for subj_gn in d.subj_gn.unique():
        for strand in ['+', '-']:
            fp = f"{args.prefix}.unambiguous.genome.{subj_gn}_{strand}.wig"
            print(f"writing to file {fp} ...")
            with open(fp, 'w') as f:
                print(f"variableStep  chrom={subj_gn}", file=f)
                for site, row in sites.loc[subj_gn, strand].iterrows():
                    print(int(site)+1, row.one_div_sites, file=f)
    # genome by subj and strands combined
    sites = ds[0].loc[d.amb == False].groupby(["subj_gn", "site_gn"])[['one_div_sites']].sum()
    for d in ds[1:]:
        sites = sites.add(d.loc[d.amb == False].groupby(["subj_gn", "site_gn"])[['one_div_sites']].sum(),
                          fill_value=0.)
    for subj_gn in d.subj_gn.unique():
        fp = f"{args.prefix}.unambiguous.genome.{subj_gn}.wig"
        print(f"writing to file {fp} ...")
        with open(fp, 'w') as f:
            print(f"variableStep  chrom={subj_gn}", file=f)
            for site, row in sites.loc[subj_gn].iterrows():
                print(int(site)+1, row.one_div_sites, file=f)
    # adapter by subj and strand
    sites = ds[0].loc[d.amb == False].groupby(["subj_ad", "strand_ad", "site_ad"])[['one_div_sites']].sum()
    for d in ds[1:]:
        sites = sites.add(d.loc[d.amb == False].groupby(["subj_ad", "strand_ad", "site_ad"])[['one_div_sites']].sum(),
                          fill_value=0.)
    for subj_ad in d.subj_ad.unique():
        for strand in ['+', '-']:
            fp = f"{args.prefix}.unambiguous.adapter.{subj_ad}_{strand}.wig"
            print(f"writing to file {fp} ...")
            with open(fp, 'w') as f:
                print(f"variableStep  chrom={subj_ad}", file=f)
                for site, row in sites.loc[subj_ad, strand].iterrows():
                    print(int(site)+1, row.one_div_sites, file=f)
    # adapter by subj and strands combined
    sites = ds[0].loc[d.amb == False].groupby(["subj_ad", "site_ad"])[['one_div_sites']].sum()
    for d in ds[1:]:
        sites = sites.add(d.loc[d.amb == False].groupby(["subj_ad", "site_ad"])[['one_div_sites']].sum(),
                          fill_value=0.)
    for subj_ad in d.subj_ad.unique():
        fp = f"{args.prefix}.unambiguous.adapter.{subj_ad}.wig"
        print(f"writing to file {fp} ...")
        with open(fp, 'w') as f:
            print(f"variableStep  chrom={subj_ad}", file=f)
            for site, row in sites.loc[subj_ad].iterrows():
                print(int(site)+1, row.one_div_sites, file=f)


if __name__ == "__main__":
    args = parse_args()
    main()