#!/usr/bin/env python

import argparse
import sys

def merge(vcf, cadd, header_only, print_header):
    for cadd_line in cadd:
        if cadd_line.startswith('#'):
            continue
        cadd_fields = cadd_line.rstrip().split('\t')

        for vcf_line in vcf:
            if vcf_line.startswith('##'):
                if header_only or print_header:
                    sys.stdout.write(vcf_line)
                continue
            elif vcf_line.startswith('#'):
                if header_only or print_header:
                    sys.stdout.write('##INFO=<ID=CADD,Number=A,Type=Float,Description="CADD score">\n')
                    sys.stdout.write('##INFO=<ID=CADD_RAW,Number=A,Type=Float,Description="Raw CADD score">\n')
                # We will always output the actual header line so we can paste in the samples
                if header_only:
                    sys.exit(0)
                sys.stdout.write(vcf_line)
            else:
                vcf_fields = vcf_line.rstrip().split('\t', 8)
                if (vcf_fields[4].startswith('<') and vcf_fields[4].endswith('>')) or vcf_fields[4]=='*':
                    vcf_fields[7] = ';'.join([vcf_fields[7], 'CADD=.', 'CADD_RAW=.'])
                    sys.stdout.write('\t'.join(vcf_fields))
                    sys.stdout.write('\n')
                elif (vcf_fields[0] == cadd_fields[0]
                        and vcf_fields[1] == cadd_fields[1]
                        and vcf_fields[3] == cadd_fields[2]
                        and vcf_fields[4] == cadd_fields[3]):
                    vcf_fields[7] = ';'.join([vcf_fields[7], 'CADD={0}'.format(cadd_fields[5]), 'CADD_RAW={0}'.format(cadd_fields[4])])
                    sys.stdout.write('\t'.join(vcf_fields))
                    sys.stdout.write('\n')
                    break
                else:
                    sys.stderr.write('HORRIBLE THINGS HAVE HAPPENED!!\n')
                    sys.exit(1)
    for vcf_line in vcf:
        vcf_fields = vcf_line.rstrip().split('\t', 8)
        if (vcf_fields[4].startswith('<') and vcf_fields[4].endswith('>')) or vcf_fields[4]=='*':
            vcf_fields[7] = ';'.join([vcf_fields[7], 'CADD=.', 'CADD_RAW=.'])
            sys.stdout.write('\t'.join(vcf_fields))
            sys.stdout.write('\n')
        else:
            sys.stderr.write('HORRIBLE THINGS HAVE HAPPENED!!\n')
            sys.exit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='vcf_file', nargs='?', help='VCF')
    parser.add_argument(dest='cadd_file', nargs='?', help='CADD')
    parser.add_argument('-H','--header-only', action='store_true', default=False, dest='header_only', help='print new header and exit')
    parser.add_argument('--no-header', action='store_false', default=True, dest='print_header', help="don't print the header")

    args = parser.parse_args()

    with open(args.vcf_file) as vcf, open(args.cadd_file) as cadd:
        merge(vcf, cadd, args.header_only, args.print_header)
