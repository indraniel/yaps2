#!/usr/bin/env python

from __future__ import print_function, division
import sys, os, datetime

if 'VIRTUAL_ENV' in os.environ:
    print('found a virtualenv -- activating: {}'.format(os.environ['VIRTUAL_ENV']), file=sys.stderr)
    activation_script = os.path.join(os.environ['VIRTUAL_ENV'], 'bin', 'activate_this.py')
    execfile(activation_script, dict(__file__=activation_script))

import click
from cyvcf2 import VCF, Writer

def log(msg):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %T")
    print('[{}] {}'.format(timestamp, msg), file=sys.stderr)

def is_float(s):
    return True if isinstance(s, float) else False

def is_int(s):
    return True if isinstance(s, (int, long) ) else False

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    rev_comp = unicode("".join(complement.get(base, base) for base in reversed(seq)))
    return rev_comp

def annotation_type_headers(annotation_type):
    headers = {
        'cadd' : [
            {
                'ID' : 'CADD',
                'Number' : 'A',
                'Type' : 'Float',
                'Description' : 'CADD score'
            },
            {
                'ID' : 'CADD_RAW',
                'Number' : 'A',
                'Type' : 'Float',
                'Description' : 'Raw CADD score'
            }
        ],
    }

    return headers[annotation_type]

def annotation_type_info_fields(annotation_type):
    fields = {
        'cadd' : ['CADD', 'CADD_RAW'],
    }

    return fields[annotation_type]

def create_b37_annotation_dictionary(b37_vcf, annotation_type):
    data = {}
    fields = annotation_type_info_fields(annotation_type)

    vcf = VCF(b37_vcf)
    for variant in vcf:
	field_data = {}
        ref = variant.REF
        alt = ','.join(variant.ALT)
        chrom = variant.INFO.get('OriginalContig')
        pos_b38 = variant.INFO.get('OriginalStart')
	for f in fields:
            value = variant.INFO.get(f, '.')
            field_data[f] = value if value else '.'

        if (chrom is None) or (pos_b38 is None):
            msg = ("In '{}' found no 'OriginalContig' and/or 'OriginalStart' "
                   "INFO attribute on "
                   "CHROM: '{}' "
                   "| POS: '{}' "
                   "| REF: '{}' "
                   "| ALT: '{}'")
            msg.format(b37_vcf, variant.CHROM, variant.POS, ref, alt)
            raise RuntimeError(msg)

        key = (chrom, pos_b38, ref, alt)
        data[key] = field_data

    return data

def update_annotations(variant, new_annotations):
    for field in new_annotations:
        fmt = '{:.4f}' if is_float(new_annotations[field]) else '{}'
        variant.INFO[field] = fmt.format(new_annotations[field])
    return variant

def update_variant(variant, b37_annotations, anno_fields, autofill):
    chrom = unicode(variant.CHROM)
    pos = unicode(variant.POS)
    ref = unicode(variant.REF)
    alt = unicode(','.join(variant.ALT))

    key = (chrom, pos, ref, alt)
    if key in b37_annotations:
        variant = update_annotations(variant, b37_annotations[key])
        return variant

    # try alternative keys based on the reverse complements
    revcomp_ref = reverse_complement(ref)
    revcomp_alt = ','.join([reverse_complement(x) for x in variant.ALT])
    key = (chrom, pos, revcomp_ref, revcomp_alt)
    if key in b37_annotations:
        variant = update_annotations(variant, b37_annotations[key])
        return variant

    # last resort, if we need to always autofill the fields
    if autofill:
        new_annotations = { x : '.' for x in anno_fields }
        variant = update_annotations(variant, new_annotations)
        return variant

    # otherwise do nothing
    return variant

def unliftover_vcf(b38_vcf, b37_vcf, annotation_type, auto_fill):
    new_info_headers = annotation_type_headers(annotation_type)
    new_annotation_fields = annotation_type_info_fields(annotation_type)

    log("Collecting the build 37 vcf annotation information")
    b37_annotations = create_b37_annotation_dictionary(b37_vcf, annotation_type)

    log("Processing the build38 vcf")
    vcf = VCF(b38_vcf)

    for info_hdr in new_info_headers:
        vcf.add_info_to_header(info_hdr)

    out = Writer('-', vcf)

    for variant in vcf:
        variant = update_variant(variant, b37_annotations, new_annotation_fields, auto_fill)
        out.write_record(variant)

    out.close()
    log("All Done!")

@click.command()
@click.option('--b38-vcf', required=True, type=click.Path(exists=True),
        help="the original b38-based vcf to integrate annotations to")
@click.option('--annotated-b37-vcf', required=True, type=click.Path(exists=True),
        help=("the b37-based vcf containing the desired annotations "
              "& OriginalContig/OriginalStart INFO fields"))
@click.option('--annotation-type',
              required=True,
              type=click.Choice(['cadd']),
              help="the type of annotation being incorporated")
@click.option('--auto-fill', is_flag=True,
              help="ensure the annotation are always populated. Insert 'FIELD=.' if empty")
def main(b38_vcf, annotated_b37_vcf, annotation_type, auto_fill):
    try:
        unliftover_vcf(b38_vcf, annotated_b37_vcf, annotation_type, auto_fill)
        return 0
    except Exception, err:
        log('[err]: {}'.format(err))
        return 1

if __name__ == "__main__":
    main()
