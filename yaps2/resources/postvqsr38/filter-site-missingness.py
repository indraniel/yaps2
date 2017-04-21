#!/usr/bin/env python

from __future__ import print_function, division
import itertools, re, sys, os

if 'VIRTUAL_ENV' in os.environ:
    print('found a virtualenv -- activating: {}'.format(os.environ['VIRTUAL_ENV']))
    activation_script = os.path.join(os.environ['VIRTUAL_ENV'], 'bin', 'activate_this.py')
    execfile(activation_script, dict(__file__=activation_script))

import click
from cyvcf2 import VCF, Writer
from pandas import Series, DataFrame
import pandas as pd

def variant_missing_criteria(threshold, variant_pct):
    if variant_pct < threshold:
        return 'pass'
    else:
        return 'fail'

def compute_missingness(variant):
    alleles = Series(list(itertools.chain.from_iterable([re.split('[/|]', x) for x in variant.gt_bases])))
    (missing, total) = ( alleles[alleles == '.'].size, alleles.size )
    missingness = (missing/total) * 100
    return (missingness, missing, total)

def update_variant(variant, verdict, missing_pct):
    if verdict == "pass":
        if not variant.FILTER:
            variant.FILTER = ['PASS']
    else:
        if (not variant.FILTER) or (variant.FILTER == 'PASS'):
            variant.FILTER = ['MISSING']
        else:
            variant.FILTER = [ variant.FILTER, 'MISSING' ]

    variant.INFO['MISSINGPCT'] = '{:.3f}'.format(missing_pct)

    return variant

def mark_missing_sites(vcffile, region, missing_threshold, soft_filter):
    vcf = VCF(vcffile)
    header_param_id = {
        'ID' : 'MISSING',
        'Description' : 'failed variant site missingness threshold ({} %)'.format(missing_threshold)
    }
    header_param_info = {
        'ID' : 'MISSINGPCT',
        'Description' : 'site missingness percentage',
        'Type' : 'Float',
        'Number' : '1'
    }
    vcf.add_filter_to_header(header_param_id)
    vcf.add_info_to_header(header_param_info)
    out = Writer('-', vcf)
    (total_sites, noted_sites) = (0, 0)

    for variant in vcf(region):
        total_sites += 1
        (missing_pct, missing, total) = compute_missingness(variant)
        verdict = variant_missing_criteria(missing_threshold, missing_pct)
        variant = update_variant(variant, verdict, missing_pct)
        if verdict == "pass":
            noted_sites += 1
            out.write_record(variant)
        elif verdict == "fail" and soft_filter:
            out.write_record(variant)

    out.close()
    msg = "After filtering, passed {} out of a possible {} Sites ({})"
    msg = msg.format(noted_sites, total_sites, 'pass')
    print(msg, file=sys.stderr)

@click.command()
@click.option('--missing-threshold', type=float, default=2.0,
        help="the missingness threshold to discard [default: 2.0]")
@click.option('--soft', is_flag=True, default=False,
        help="soft filtering -- keep all variants and update FILTER field [default: False]")
@click.option('--region', default=None, type=click.STRING,
        help="a chromosome region to limit to [default: None]")
@click.argument('vcfs', nargs=-1, type=click.Path())
def main(missing_threshold, region, soft, vcfs):
    for vcf in vcfs:
        mark_missing_sites(vcf, region, missing_threshold, soft)
    print("All Done!", file=sys.stderr)

if __name__ == "__main__":
    main()
