#!/usr/bin/env python

from __future__ import print_function, division
import sys, os, datetime, gzip

if 'VIRTUAL_ENV' in os.environ:
    print('found a virtualenv -- activating: {}'.format(os.environ['VIRTUAL_ENV']), file=sys.stderr)
    activation_script = os.path.join(os.environ['VIRTUAL_ENV'], 'bin', 'activate_this.py')
    execfile(activation_script, dict(__file__=activation_script))

import click
from cyvcf2 import VCF, Writer

def log(msg):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %T")
    print('[-- {} --] {}'.format(timestamp, msg), file=sys.stderr)

def is_float(s):
    return True if isinstance(s, float) else False

def annotation_info_headers():
    headers = [
        {
            'ID' : 'CADD',
            'Number' : 'A',
            'Type' : 'Float',
            'Description' : 'CADD score',
        },
        {
            'ID' : 'CADD_RAW',
            'Number' : 'A',
            'Type' : 'Float',
            'Description' : 'Raw CADD score',
        }
    ]
    return headers

def create_CADD_annotation_dictionary(tsv):
    data = {}

    with gzip.GzipFile(tsv, 'rb') as f:
        for line in f:
            if line.startswith('#'):
                continue
            (chrom, pos, ref, alt, raw_score, phred_score) = line.rstrip().split("\t")
            chrom = unicode(chrom)
            pos = unicode(pos)
            ref = unicode(ref)
            alt = unicode(alt)
            key = (chrom, pos, ref, alt)
            data[key] = { 'CADD' : phred_score, 'CADD_RAW' : raw_score }

    return data

def update_annotations(variant, cadd_score, raw_cadd_score):
    cadd_score_fmt = '{:.3f}' if is_float(cadd_score) else '{}'
    raw_cadd_score_fmt = '{:.6f}' if is_float(raw_cadd_score) else '{}'

    variant.INFO['CADD'] = cadd_score_fmt.format(cadd_score)
    variant.INFO['CADD_RAW'] = raw_cadd_score_fmt.format(raw_cadd_score)

    return variant

def update_variant(variant, cadd_annotations):
    chrom = unicode(variant.CHROM)
    pos = unicode(variant.POS)
    ref = unicode(variant.REF)
    alt = unicode(','.join(variant.ALT))

    key = (chrom, pos, ref, alt)
    if key in cadd_annotations:
        cadd_score = cadd_annotations[key]['CADD']
        raw_cadd_score = cadd_annotations[key]['CADD_RAW']
    else:
        cadd_score = '.'
        raw_cadd_score = '.'

    update_annotations(variant, cadd_score, raw_cadd_score)
    return variant

def merge(in_vcf, cadd_tsv):
    new_headers = annotation_info_headers()

    log("Collecting the CADD annotation information")
    cadd_annotations = create_CADD_annotation_dictionary(cadd_tsv)

    log("Processing the build38 vcf")
    vcf = VCF(in_vcf)

    for info_hdr in new_headers:
        vcf.add_info_to_header(info_hdr)

    out = Writer('-', vcf)

    for variant in vcf:
        variant = update_variant(variant, cadd_annotations)
        out.write_record(variant)

    out.close()
    log("All Done!")

@click.command()
@click.option('--in-vcf', required=True, type=click.Path(exists=True),
        help="the original input vcf to integrate annotations into")
@click.option('--cadd-tsv', required=True, type=click.Path(exists=True),
        help="the gzipped tsv file produced by CADD's score.sh")
def main(in_vcf, cadd_tsv):
    merge(in_vcf, cadd_tsv)

if __name__ == "__main__":
    main()
