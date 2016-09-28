#!/usr/bin/env python

from __future__ import print_function, division
import os

if 'VIRTUAL_ENV' in os.environ:
    print('found a virtualenv -- activating: {}'.format(os.environ['VIRTUAL_ENV']))
    activation_script = os.path.join(os.environ['VIRTUAL_ENV'], 'bin', 'activate_this.py')
    execfile(activation_script, dict(__file__=activation_script))

import click

def aggregate(basedir, partitions, chroms):
    stats = {}
    families = {}
    for bucket in partitions:
        stats[bucket] = {}
        for chrom in chroms:
            inspect_dir = os.path.join(basedir, str(bucket), str(chrom))
            data = process_stats_bucket_chrom_dir(inspect_dir)
            for family in data:
                families.setdefault(family, 1)
                prior_counts = stats[bucket].setdefault(family, {'variants' : 0, 'mie' : 0 })
                prior_counts['variants'] = prior_counts['variants'] + data[family]['variants']
                prior_counts['mie'] = prior_counts['mie'] + data[family]['mie']
                stats[bucket][family] = prior_counts

    # ensure all the families are in every bucket
    for bucket in stats:
        for family in families:
            if family not in stats[bucket]:
                stats[bucket].setdefault(family, {'variants' : 0, 'mie' : 0 })

    return stats

def print_stats_table(stats, category, method, out_file):
    header = ['category', 'method', 'bucket', 'family', 'variants', 'mie']
    with open(out_file, 'w') as f:
        print("\t".join(header), file=f)
        for bucket in sorted(stats.keys()):
            for family in sorted(stats[bucket].keys()):
                line = [
                    category,
                    method,
                    bucket,
                    family,
                    stats[bucket][family]['variants'],
                    stats[bucket][family]['mie'],
                ]
                line = [ str(i) for i in line ]
                print("\t".join(line), file=f)

def process_stats_bucket_chrom_dir(bucket_dir):
    print("processing: {}".format(bucket_dir))
    empty_vcf = os.path.join(bucket_dir, 'empty-vcf')
    var_txt = os.path.join(bucket_dir, 'unfiltered.mie.var.txt')

    if os.path.isfile(empty_vcf):
        return {}
    elif os.path.isfile(var_txt):
        stats = process_unfiltered_mie_var_txt_file(var_txt)
        return stats
    else:
        msg = ("Summary MIE file (unfiltered.mie.var.txt) or "
               "'empy-vcf' file not found!")
        raise RuntimeError(msg)

def process_unfiltered_mie_var_txt_file(var_txt):
    stats = {}
    with open(var_txt, 'r') as f:
        for line in f:
            # skip header
            if '\t' in line:
                continue
            line = line.rstrip()
            (family, variants, mie) = line.split(' ')
            stats[family] = { 'variants' : int(variants), 'mie' : int(mie) }
    return stats

def calculate_mie_stats(input_dir, category, method, tranches, percentiles, chroms):
    basedir = os.path.join(input_dir, category, method)
    if method == 'tranche':
        stats = aggregate(basedir, tranches, chroms)
    else:
        stats = aggregate(basedir, percentiles, chroms)
    return stats

def ensure_output_dir(out_file):
    basedir = os.path.dirname(out_file)
    if not os.path.isdir(basedir):
        os.makedirs(basedir)

def cli_setup_tranches(ctx, param, value):
    if value is None:
        return tuple(range(1,4))
    else:
        return tuple(int(i.strip()) for i in value.split(','))

def cli_setup_percentiles(ctx, param, value):
    if value is None:
        return tuple(range(10, 110, 10))
    else:
        return tuple(int(i.strip()) for i in value.split(','))

def cli_setup_chroms(ctx, param, value):
    if value is None:
        return tuple(range(1,23))
    else:
        return tuple(int(i.strip()) for i in value.split(','))

@click.command()
@click.option('--input-dir', required=True, type=click.Path(exists=True),
              help='The base input directory containing the relevant MIE stats')
@click.option('--output-file', required=True, type=click.Path(),
              help='The base input directory containing the relevant MIE stats')
@click.option('--category', required=True, type=click.Choice(['snps', 'indels']),
              help='The category of stats')
@click.option('--method', required=True, type=click.Choice(['tranche', 'percentile']),
              help='The category of stat breakdowns')
@click.option('--tranches', default=None, type=click.STRING, callback=cli_setup_tranches,
              help='A comma delimited set of tranches to inspect [default="1,2,3"]')
@click.option('--percentiles', default=None, type=click.STRING, callback=cli_setup_percentiles,
              help=('A comma delimited set of percentiles to inspect '
                    '[default="10,20,30,40,50,60,70,80,90,100"]'))
@click.option('--chroms', default=None, type=click.STRING, callback=cli_setup_chroms,
              help=('A comma delimited set of chromosomes to inspect '
                    '[default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"]'))
def main(input_dir, output_file, category, method, tranches, percentiles, chroms):
    stats = calculate_mie_stats(input_dir, category, method, tranches, percentiles, chroms)
    ensure_output_dir(output_file)
    print_stats_table(stats, category, method, output_file)

if __name__ == '__main__':
    main()
