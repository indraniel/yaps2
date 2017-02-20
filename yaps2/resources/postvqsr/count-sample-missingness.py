#!/usr/bin/env python

from __future__ import print_function, division
import sys, os, json, gzip

if 'VIRTUAL_ENV' in os.environ:
    print('found a virtualenv -- activating: {}'.format(os.environ['VIRTUAL_ENV']))
    activation_script = os.path.join(os.environ['VIRTUAL_ENV'], 'bin', 'activate_this.py')
    execfile(activation_script, dict(__file__=activation_script))

import click
from cyvcf2 import VCF, Writer
import numpy as np

def calculate_sample_missingness(vcffile):
    vcf = VCF(vcffile)
    missing_counts = np.zeros(len(vcf.samples)).astype(np.uint64)
    total_passing_variants = 0

    for variant in vcf:
        if (variant.FILTER is not None) and (variant.FILTER != 'PASS'): continue
        total_passing_variants += 1
        genotypes = variant.gt_bases
        mask = (genotypes == './.') | (genotypes == '.|.')
        missing_counts = missing_counts + mask.astype(np.uint64)

    sample_missingness = dict(zip(vcf.samples, missing_counts))
    stats = {
        'total_pass_variants' : total_passing_variants,
        'missingness_counts' : sample_missingness,
    }
    return stats

def merge_stats(vcfstats, totals):
    merged = { 'total_pass_variants' : 0, 'missingness_counts' : {} }
    merged['total_pass_variants'] = vcfstats['total_pass_variants'] + totals.get('total_pass_variants', 0)

    if 'missingness_counts' in totals:
        vcfc = vcfstats['missingness_counts']
        totalc = totals['missingness_counts']
        mergedc = merged['missingness_counts']

        for sample in vcfc:
            mergedc[sample] = vcfc[sample] + totalc.get(sample, 0)

        # for now assume all the input vcfs have the exact same samples
        # error out if there are samples in the overall total dictionary
        # that are not in the individual vcf stats
        error_samples = [sample for sample in totalc if sample not in vcfc]
        if len(error_samples > 0):
            msg = ( '[err] The following {} samples are unaccounted for. '
                    'Please investigate!\n{}\n' )
            sys.exit(msg.format(len(error_samples), error_samples))
    else:
        merged['missingness_counts'] = vcfstats['missingness_counts']

    return merged

def dump_stats(outfile, stats):
    with open(outfile, 'w') as f:
        data = json.dumps(eval(str(stats)), sort_keys=True, indent=4, separators=(',', ': '))
        print(data, file=f)

@click.command()
@click.option('--out', default="sample-missingness.out", type=click.Path(),
        help="an output file to write site stats to [default: 'sample-missingness.out']")
@click.argument('vcfs', nargs=-1, type=click.Path())
def main(out, vcfs):
    totals = {}
    for vcf in vcfs:
        stats = calculate_sample_missingness(vcf)
        totals = merge_stats(stats, totals)
    dump_stats(out, totals)
    print("All Done!", file=sys.stderr)

if __name__ == "__main__":
    main()
