#!/usr/bin/env python

from __future__ import print_function, division
import sys, os, json, datetime, time

if 'VIRTUAL_ENV' in os.environ:
    print('found a virtualenv -- activating: {}'.format(os.environ['VIRTUAL_ENV']))
    activation_script = os.path.join(os.environ['VIRTUAL_ENV'], 'bin', 'activate_this.py')
    execfile(activation_script, dict(__file__=activation_script))

import click

def logit(msg):
    ts = time.strftime("[ %Y-%m-%d %T ]", datetime.datetime.now().timetuple())
    fullmsg = "{} {}".format(ts, msg)
    print(fullmsg)
    sys.stdout.flush()
    sys.stderr.flush()

def get_data(jsonfile):
    with open(jsonfile) as f:
        data = json.load(f)
    return data

def merge_stats(jsonfile, totals):
    stats = get_data(jsonfile)

    if stats['total_pass_variants'] == 0:
        return totals

    merged = { 'total_pass_variants' : 0, 'missingness_counts' : {} }
    merged['total_pass_variants'] = stats['total_pass_variants'] + totals.get('total_pass_variants', 0)

    if 'missingness_counts' in totals:
        statsc = stats['missingness_counts']
        totalc = totals['missingness_counts']
        mergedc = merged['missingness_counts']

        for sample in statsc:
            # for now assume all the input json files have the exact same samples
            # error out if a sample is not found in the overall total dictionary
            if sample not in totalc:
                msg = ( "[err] Sample '{}' in '{}' is unaccounted for "
                        "in other json documents. "
                        "Please investigate!" )
                sys.exit(msg.format(sample, jsonfile))
            mergedc[sample] = statsc[sample] + totalc[sample]
        merged['missingness_counts'] = mergedc
    else:
        merged['missingness_counts'] = stats['missingness_counts']

    return merged

def dump_stats(outfile, data):
    total = data['total_pass_variants']
    counts = data['missingness_counts']

    with open(outfile, 'w') as f:
        headers = ('#SAMPLE', 'MISSING', 'TOTAL_PASSING_VARIANTS', 'PCT_MISSING')
        print("{:30} {:15} {:30} {}".format(*headers), file=f)

        fmt = "{:<30} {:<15} {:<30} {:<.4f}"
        for sample in sorted(counts.keys()):
            missing = counts[sample]
            pct_missing = ( missing / total ) * 100.0
            print(fmt.format(sample, missing, total, pct_missing), file=f)

@click.command()
@click.option('--out', default="sample-missingness.out", type=click.Path(),
        help="an output file to write site stats to [default: 'sample-missingness.out']")
@click.argument('jsondocs', nargs=-1, type=click.Path())
def main(out, jsondocs):
    totals = {}
    for jsonfile in jsondocs:
        logit("Merging {}".format(jsonfile))
        totals = merge_stats(jsonfile, totals)
    logit("Dumping overall sample missingness statistics to {}".format(out))
    dump_stats(out, totals)
    logit("All Done!")

if __name__ == "__main__":
    main()
