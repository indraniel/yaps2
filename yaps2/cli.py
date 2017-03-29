from __future__ import print_function, division
import signal, sys, os

import click

from .version import __version__

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=__version__)
def cli():
    '''A collection of CCDG related VCF data processing and QC analysis pipelines.'''
    # to make this script/module behave nicely with unix pipes
    # http://newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

@cli.command(short_help="post-VQSR data pipeline")
@click.option('--workspace', required=True, type=click.Path(),
              help='A directory to place outputs into')
@click.option('--job-db', default=None, type=click.Path(),
              help="Path to LSF job sqlite DB [default='<workspace>/.job_queue.db']")
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
@click.option('--docker/--no-docker', default=False,
              help='Use the "docker-ize" pipeline [default=False or --no-docker]')
@click.option('--skip-confirm', default=False, is_flag=True,
              help='Do not prompt when resuming or restarting a pipeline [default=False]')
@click.option('--task-flush', default=False, is_flag=True,
              help='Update the task database table as soon as a job is submitted [default=False]')
def postvqsr(job_db, input_vcfs, project_name, email, workspace, drm, restart, docker, skip_confirm, task_flush):
    from yaps2.pipelines.postvqsr import Config, Pipeline
    config = Config(job_db, input_vcfs, project_name, email, workspace, docker)
    workflow = Pipeline(config, drm, restart, skip_confirm)
    workflow.run(task_flush)

@cli.command(short_help="Mendelian Inheritance Error [MIE] Analysis Pipeline")
@click.option('--workspace', required=True, type=click.Path(),
              help='A directory to place outputs into')
@click.option('--job-db', default=None, type=click.Path(),
              help="Path to LSF job sqlite DB [default='<workspace>/.job_queue.db']")
@click.option('--input-vcfs', required=True, type=click.Path(exists=True),
              help='A file of chromosomal VCFs to process')
@click.option('--samples', required=True, type=click.Path(exists=True),
              help='A list of samples to restrict to [sample-name delimited by newlines]')
@click.option('--plink-fam', required=True, type=click.Path(exists=True),
              help='A plink .fam file describing the sample relationships')
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
def mie(job_db, input_vcfs, percentiles, samples, tranches, plink_fam,
        project_name, email, workspace, drm, restart):
    from yaps2.pipelines.mie import Config, Pipeline
    config = Config(
        job_db, input_vcfs,
        percentiles, samples, tranches, plink_fam,
        project_name, email, workspace
    )
    workflow = Pipeline(config, drm, restart)
    workflow.run()

@cli.command(short_help="Principal Component Analysis [PCA] Pipeline")
@click.option('--workspace', required=True, type=click.Path(),
              help='A directory to place outputs into')
@click.option('--job-db', default=None, type=click.Path(),
              help="Path to LSF job sqlite DB [default='<workspace>/.job_queue.db']")
@click.option('--input-vcfs', required=True, type=click.Path(exists=True),
              help='A file of chromosomal VCFs to process')
@click.option('--project-name', default='yaps2.default', type=click.STRING,
              help='A prefix used to name batch jobs')
@click.option('--vqslod-threshold', required=True, type=click.FLOAT,
              help='A minimum VQSLOD tranche threshold to choose SNPs from [min. snp tranche 2 VQSLOD point]')
@click.option('--email', default=None, type=click.STRING,
              help='An email used to notify about batch jobs [default=userid@genome.wustl.edu]')
@click.option('--drm', default='lsf', type=click.Choice(['local', 'lsf']),
              help='Job Mode -- [default=lsf]')
@click.option('--restart/--no-restart', default=False,
              help='Restart Pipeline from scratch')
def pca(job_db, input_vcfs, project_name, email, workspace, vqslod_threshold, drm, restart):
    from yaps2.pipelines.pca import Config, Pipeline
    config = Config(job_db, input_vcfs, project_name, email, workspace, vqslod_threshold)
    workflow = Pipeline(config, drm, restart)
    workflow.run()

@cli.command(name='b38-realign', short_help="Re-align raw sequence data with Build 38 and speedseq")
@click.option('--workspace', required=True, type=click.Path(),
              help='A directory to place outputs into')
@click.option('--job-db', default=None, type=click.Path(),
              help="Path to LSF job sqlite DB [default='<workspace>/.job_queue.db']")
@click.option('--project-name', default='yaps2.default', type=click.STRING,
              help='A prefix used to name batch jobs')
@click.option('--input-sample-bams', required=True, type=click.Path(),
              help='A JSON file containing the sample/bam data to be processed')
@click.option('--email', default=None, type=click.STRING,
              help='An email used to notify about batch jobs [default=userid@genome.wustl.edu]')
@click.option('--drm', default='lsf', type=click.Choice(['local', 'lsf']),
              help='Job Mode -- [default=lsf]')
@click.option('--drm-job-group', required=True, type=click.STRING,
              help='An LSF job group to control cluster usage')
@click.option('--restart/--no-restart', default=False,
              help='Restart Pipeline from scratch')
def b38_realign(job_db, input_sample_bams, project_name, email, workspace, drm, drm_job_group, restart):
    from yaps2.pipelines.b38 import Config, Pipeline
    config = Config(job_db, input_sample_bams, project_name, email, workspace, drm_job_group)
    workflow = Pipeline(config, drm, restart)
    workflow.run()
