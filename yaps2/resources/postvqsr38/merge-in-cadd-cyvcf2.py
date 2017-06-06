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
    return variant, key

def ensure_cadd_completed_successfully(in_vcf_file, cadd_tsv_file, vcf_set, cadd_set):
    # In theory, with AC=0, symbolic deletions, and the unplaced contigs all removed,
    # the number of (chrom, pos, ref, alt) combinations should be the in the input vcf
    # and output CADD tsv file
    num_in_vcf_variants = len(vcf_set)
    num_in_cadd_tsv_variants = len(cadd_set)
    if num_in_vcf_variants != num_in_cadd_tsv_variants:
        diff = vcf_set.symmetric_difference(cadd_set)
        msg = ("[err] Found {} variants that are "
               "not common between CADD's input and outputs")
        log( msg.format(len(diff)) )
        log( "Input CADD VCF  : '{}'".format(in_vcf_file) )
        log( "Output CADD TSV : '{}'".format(cadd_tsv_file) )
        log( "The troublesome variants are:" )

        detail = '    CHROM: {} | POS: {} | REF: {} | ALT: {}'
        for variant in diff:
            (chrom, pos, ref, alt) = variant
            print(detail.format(chrom, pos, ref, alt), file=sys.stderr)

        raise RuntimeError( msg.format(len(diff)) )

    log("Successfully passed check")

def merge(in_vcf, cadd_tsv):
    new_headers = annotation_info_headers()

    log("Collecting the CADD annotation information")
    cadd_annotations = create_CADD_annotation_dictionary(cadd_tsv)

    log("Processing the build37 vcf")
    vcf = VCF(in_vcf)

    for info_hdr in new_headers:
        vcf.add_info_to_header(info_hdr)

    out = Writer('-', vcf)

    in_vcf_variants = set()
    for variant in vcf:
        (variant, key) = update_variant(variant, cadd_annotations)
        in_vcf_variants.add(key)
        out.write_record(variant)

    out.close()

    log("Checking whether CADD completed correctly")
    ensure_cadd_completed_successfully(in_vcf, cadd_tsv, in_vcf_variants, frozenset(list(cadd_annotations.keys())) )

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
