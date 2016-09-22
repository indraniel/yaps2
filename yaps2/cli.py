from __future__ import print_function, division
import signal, sys, os

import click

from .version import __version__

@click.group()
@click.version_option(version=__version__)
def cli():
    # to make this script/module behave nicely with unix pipes
    # http://newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

@cli.command()
@click.option('--workspace', required=True, type=click.Path(),
              help='A directory to place outputs into')
@click.option('--job-db', default=None, type=click.Path(),
              help="Path to LSF job sqlite DB [default='<workspace>/job_queue.db']")
@click.option('--input-vcfs', required=True, type=click.Path(exists=True),
              help='A file of chromosomal VCFs to process')
@click.option('--project-name', default='yaps2.default', type=click.STRING,
              help='A prefix used to name batch jobs')
@click.option('--email', default=None, type=click.STRING,
              help='An email used to notify about batch jobs [default=userid@genome.wustl.edu]')
@click.option('--drm', default='lsf', type=click.Choice(['local', 'lsf']),
              help='Job Mode -- [default=lsf]')
@click.option('--restart/--no-restart', default=False,
              help='Restart Pipeline from scratch')
def postvqsr(job_db, input_vcfs, project_name, email, workspace, drm, restart):
    from yaps2.pipelines.postvqsr import Config, Pipeline
    config = Config(job_db, input_vcfs, project_name, email, workspace)
    workflow = Pipeline(config, drm, restart)
    workflow.run()

@cli.command()
@click.option('--workspace', required=True, type=click.Path(),
              help='A directory to place outputs into')
@click.option('--job-db', default=None, type=click.Path(),
              help="Path to LSF job sqlite DB [default='<workspace>/job_queue.db']")
@click.option('--input-vcfs', required=True, type=click.Path(exists=True),
              help='A file of chromosomal VCFs to process')
@click.option('--samples', required=True, type=click.Path(exists=True),
              help='A list of samples to restrict to [sample-name delimited by newlines]')
@click.option('--percentiles', required=True, type=click.Path(exists=True),
              help='A tsv of category/percentile/min-vqslod/max-vqslod lines')
@click.option('--tranches', required=True, type=click.Path(exists=True),
              help='A tsv of category/tranche/min-vqslod lines')
@click.option('--project-name', default='yaps2.default', type=click.STRING,
              help='A prefix used to name batch jobs')
@click.option('--email', default=None, type=click.STRING,
              help='An email used to notify about batch jobs [default=userid@genome.wustl.edu]')
@click.option('--drm', default='lsf', type=click.Choice(['local', 'lsf']),
              help='Job Mode -- [default=lsf]')
@click.option('--restart/--no-restart', default=False,
              help='Restart Pipeline from scratch')
def mie(job_db, input_vcfs, percentiles, samples, tranches,
        project_name, email, workspace, drm, restart):
    from yaps2.pipelines.mie import Config, Pipeline
    config = Config(
        job_db, input_vcfs,
        percentiles, samples, tranches,
        project_name, email, workspace
    )
    workflow = Pipeline(config, drm, restart)
    workflow.run()
