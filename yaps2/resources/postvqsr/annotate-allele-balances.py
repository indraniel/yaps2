#!/usr/bin/env python

from __future__ import print_function, division
import sys, os, datetime

if 'VIRTUAL_ENV' in os.environ:
    print('found a virtualenv -- activating: {}'.format(os.environ['VIRTUAL_ENV']), file=sys.stderr)
    activation_script = os.path.join(os.environ['VIRTUAL_ENV'], 'bin', 'activate_this.py')
    execfile(activation_script, dict(__file__=activation_script))

import click
from cyvcf2 import VCF, Writer
import numpy as np

def log(msg):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %T")
    print('[{}] {}'.format(timestamp, msg), file=sys.stderr)

def is_biallelic(variant):
    return True if len(variant.ALT) == 1 else False

def compute_allelic_balances(variant):
    # initial values
    (het_ab, het_hom_alt_ab) = (0.0 , 0.0)

    # gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
    het_mask = variant.gt_types == 1
    ref_counts = np.sum(variant.format('AD')[het_mask][:,0])
    alt_counts = np.sum(variant.format('AD')[het_mask][:,1])
    total_counts = alt_counts + ref_counts
    if total_counts != 0:
        het_ab = alt_counts / total_counts

    het_hom_alt_mask = (variant.gt_types == 1) | (variant.gt_types == 3)
    ref_counts = np.sum(variant.format('AD')[het_hom_alt_mask][:,0])
    alt_counts = np.sum(variant.format('AD')[het_hom_alt_mask][:,1])
    total_counts = alt_counts + ref_counts
    if total_counts != 0:
        het_hom_alt_ab = alt_counts / total_counts

    return (het_ab, het_hom_alt_ab)

def update_variant(variant, het_ab, het_hom_alt_ab):
    variant.INFO['HetAB'] = '{:.4f}'.format(het_ab)
    variant.INFO['HetHomAltAB'] = '{:.4f}'.format(het_hom_alt_ab)
    return variant

def annotate_allelic_balance(vcffile, region):
    vcf = VCF(vcffile)

    header_hetab_param_info = {
        'ID' : 'HetAB',
        'Description' : 'heterozygous genotype allele balance',
        'Type' : 'Float',
        'Number' : '1'
    }

    header_het_hom_alt_ab_param_info = {
        'ID' : 'HetHomAltAB',
        'Description' : 'heterozygous + homozygous ALT genotype allele balance',
        'Type' : 'Float',
        'Number' : '1'
    }

    vcf.add_info_to_header(header_hetab_param_info)
    vcf.add_info_to_header(header_het_hom_alt_ab_param_info)
    out = Writer('-', vcf)
    (total_sites, noted_sites) = (0, 0)

    for variant in vcf(region):
        total_sites += 1
        if is_biallelic(variant):
            noted_sites += 1
            (hetab, het_hom_alt_ab) = compute_allelic_balances(variant)
            variant = update_variant(variant, hetab, het_hom_alt_ab)
        out.write_record(variant)

    out.close()
    msg = "Annotated {} out of a possible {} sites"
    msg = msg.format(noted_sites, total_sites)
    log(msg)

@click.command()
@click.option('--region', default=None, type=click.STRING,
        help="a chromosome region to limit to [default: None]")
@click.argument('vcfs', nargs=-1, type=click.Path())
def main(region, vcfs):
    for vcf in vcfs:
        log('processing: {}'.format(vcf))
        annotate_allelic_balance(vcf, region)
    log("All Done!")

if __name__ == "__main__":
    main()
