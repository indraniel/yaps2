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
                'Description' : 'CADD score',
            },
            {
                'ID' : 'CADD_RAW',
                'Number' : 'A',
                'Type' : 'Float',
                'Description' : 'Raw CADD score',
            }
        ],
        '1000G' : [
            {
                'ID' : '1KG_EAS_AF',
                'Number' : 'A',
                'Type' : 'Float',
                'Description' : 'Allele frequency in the EAS populations calculated from AC and AN, in the range (0,1)',
            },
            {
                'ID' : '1KG_EUR_AF',
                'Number' : 'A',
                'Type' : 'Float',
                'Description' : 'Allele frequency in the EUR populations calculated from AC and AN, in the range (0,1)',
            },
            {
                'ID' : '1KG_AFR_AF',
                'Number' : 'A',
                'Type' : 'Float',
                'Description' : 'Allele frequency in the AFR populations calculated from AC and AN, in the range (0,1)',
            },
            {
                'ID' : '1KG_AMR_AF',
                'Number' : 'A',
                'Type' : 'Float',
                'Description' : 'Allele frequency in the AMR populations calculated from AC and AN, in the range (0,1)',
            },
            {
                'ID' : '1KG_SAS_AF',
                'Number' : 'A',
                'Type' : 'Float',
                'Description' : 'Allele frequency in the SAS populations calculated from AC and AN, in the range (0,1)',
            },
        ]
    }

    return headers[annotation_type]

def annotation_type_info_fields(annotation_type):
    fields = {
        'cadd' : ['CADD', 'CADD_RAW'],
        '1000G': ['1KG_EAS_AF', '1KG_EUR_AF', '1KG_AFR_AF', '1KG_AMR_AF', '1KG_SAS_AF'],
    }

    return fields[annotation_type]

def create_b37_annotation_dictionary(b37_vcf, annotation_type, update_id_flag):
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

        if update_id_flag and (variant.ID is not None):
            field_data['ID'] = variant.ID

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

def update_annotations(variant, anno_fields, new_annotations):
    for field in anno_fields:
        fmt = '{:.4f}' if is_float(new_annotations[field]) else '{}'
        variant.INFO[field] = fmt.format(new_annotations[field])
    return variant

def update_vcf_id(variant, new_id):
    if new_id is not None:
        variant.ID = new_id
    return variant

def update_variant(variant, b37_annotations, anno_fields, autofill, update_id):
    chrom = unicode(variant.CHROM)
    pos = unicode(variant.POS)
    ref = unicode(variant.REF)
    alt = unicode(','.join(variant.ALT))

    key = (chrom, pos, ref, alt)
    if key in b37_annotations:
        variant = update_annotations(variant, anno_fields, b37_annotations[key])
        if update_id:
            variant = update_vcf_id(variant, b37_annotations[key].get('ID', None))
        return variant

    # try alternative keys based on the reverse complements
    revcomp_ref = reverse_complement(ref)
    revcomp_alt = ','.join([reverse_complement(x) for x in variant.ALT])
    key = (chrom, pos, revcomp_ref, revcomp_alt)
    if key in b37_annotations:
        variant = update_annotations(variant, anno_fields, b37_annotations[key])
        if update_id:
            variant = update_vcf_id(variant, b37_annotations[key].get('ID', None))
        return variant

    # last resort, if we need to always autofill the fields
    if autofill:
        new_annotations = { x : '.' for x in anno_fields }
        variant = update_annotations(variant, anno_fields, new_annotations)
        if update_id:
            variant = update_vcf_id(variant, b37_annotations[key].get('ID', None))
        return variant

    # otherwise do nothing
    return variant

def unliftover_vcf(b38_vcf, b37_vcf, annotation_type, auto_fill, update_id):
    new_info_headers = annotation_type_headers(annotation_type)
    new_annotation_fields = annotation_type_info_fields(annotation_type)

    log("Collecting the build 37 vcf annotation information")
    b37_annotations = create_b37_annotation_dictionary(b37_vcf, annotation_type, update_id)

    log("Processing the build38 vcf")
    vcf = VCF(b38_vcf)

    for info_hdr in new_info_headers:
        vcf.add_info_to_header(info_hdr)

    out = Writer('-', vcf)

    for variant in vcf:
        variant = update_variant(variant, b37_annotations, new_annotation_fields, auto_fill, update_id)
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
              type=click.Choice(['cadd', '1000G']),
              help="the type of annotation being incorporated")
@click.option('--auto-fill', is_flag=True,
              help="ensure the annotation are always populated. Insert 'FIELD=.' if empty")
@click.option('--update-id', is_flag=True,
              help="update the ID Field")
def main(b38_vcf, annotated_b37_vcf, annotation_type, auto_fill, update_id):
    try:
        unliftover_vcf(b38_vcf, annotated_b37_vcf, annotation_type, auto_fill, update_id)
        return 0
    except Exception, err:
        log('[err]: {}'.format(err))
        return 1

if __name__ == "__main__":
    main()
