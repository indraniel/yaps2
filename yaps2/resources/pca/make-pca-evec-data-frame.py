#!/usr/bin/env python

from __future__ import print_function, division
import os

if 'VIRTUAL_ENV' in os.environ:
    print('found a virtualenv -- activating: {}'.format(os.environ['VIRTUAL_ENV']))
    activation_script = os.path.join(os.environ['VIRTUAL_ENV'], 'bin', 'activate_this.py')
    execfile(activation_script, dict(__file__=activation_script))

import click

def make_data_frame(src, dst):
    with open(src, 'r') as fin:
        with open(dst, 'w') as fout:
            orig_header = fin.readline().rstrip()
            hdr_components = orig_header.split()
            new_header = ['SAMPLE']
            new_header.extend( ['PC{}'.format(i+1) for (i,item) in enumerate(hdr_components[1:])] )
            new_header.append('PHENOTYPE')
            print('\t'.join(hdr_components), file=fout)
            print('\t'.join(new_header), file=fout)
            for line in fin:
                entries = line.rstrip().split()
                print('\t'.join(entries), file=fout)

@click.command()
@click.option('--src', required=True, type=click.Path(exists=True),
              help='The raw *.pca.evec file produced by eigenstrat/smartpca')
@click.option('--out', required=True, type=click.Path(),
              help='The destination output file')
def main(src, out):
    make_data_frame(src, out)
    print("Successfully created: {}".format(out))

if __name__ == '__main__':
    main()
