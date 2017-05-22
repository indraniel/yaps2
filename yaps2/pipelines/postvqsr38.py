import os, pwd, sys
import pkg_resources
from itertools import groupby
from cosmos.api import Cosmos, Dependency, default_get_submit_args
from yaps2.utils import to_json, merge_params, natural_key, empty_gzipped_vcf

class Config(object):
    def __init__(self, job_db, input_vcf_list, project_name, email, workspace, docker):
        self.email = email
        self.db = job_db
        self.project_name = project_name
        self.rootdir = workspace
        self.docker = docker

        self.ensure_rootdir()

        if self.email is None:
            self.email = self.setup_email()

        if self.db is None:
            self.db = os.path.join(
                os.path.abspath(self.rootdir),
                '.job_queue.db'
            )

        self.vcfs = self.collect_input_vcfs(input_vcf_list)
        self.chroms = self.get_ordered_chroms()

    def ensure_rootdir(self):
        if not os.path.exists(self.rootdir):
            os.makedirs(self.rootdir)

    def setup_email(self):
        user_id = pwd.getpwuid( os.getuid() ).pw_name
        return '{}@genome.wustl.edu'.format(user_id)

    def collect_input_vcfs(self, infile):
        # expecting a tsv file of <chrom>\t<path-to-vcf-file> lines
        with open(infile, 'r') as f:
            vcfs = [ tuple(line.rstrip().split("\t")) for line in f ]
        return dict(vcfs)

    def get_ordered_chroms(self):
        chroms = sorted(self.vcfs.keys(), key=natural_key)
        return chroms

class Pipeline(object):
    def __init__(self, config, drm, restart, skip_confirm):
        self.config = config

        self.cosmos = Cosmos(
            database_url='sqlite:///{}'.format(self.config.db),
            get_submit_args=default_get_submit_args,
            default_drm=drm
        )

        self.cosmos.initdb()

        primary_logfile = os.path.join(
            self.config.rootdir,
            '{}.log'.format(self.config.project_name),
        )

        self.workflow = self.cosmos.start(
            self.config.project_name,
            primary_log_path=primary_logfile,
            restart=restart,
            skip_confirm=skip_confirm,
        )

        self.setup_pipeline()

    def setup_pipeline(self):
        self.construct_pipeline()
        self.workflow.make_output_dirs()

    def run(self, task_flush):
	# put set_successful to False if you intend to add more tasks to the
	# pipeline later
        custom_log_dir = lambda task : os.path.join(self.config.rootdir, 'logs', task.stage.name, task.uid)
        self.workflow.run(set_successful=False, log_out_dir_func=custom_log_dir, db_task_flush=task_flush)

    def construct_pipeline(self):
        # 1. calculate sample missingness (counting phase)
        count_sample_missingness_tasks = self.create_count_sample_missingness_tasks(1)
        # 1.1 calculate sample missingness (merge and calculation phase)
        calculate_sample_missingness_task = self.create_calculate_sample_missingness_task(count_sample_missingness_tasks, 1.1)
        # 2. denormalize, decompose, and uniq
        dnu_tasks = self.create_decompose_normalize_unique_tasks(2)
        # 3. remove symbolic alleles
        rsa_tasks = self.create_remove_symbolic_deletion_tasks(dnu_tasks, 3)
        # 4. filter missingness
        filter_variant_missingness_tasks = self.create_filter_variant_missingness_tasks(rsa_tasks, 4)
        # 5. annotate allele balances
        allele_balance_annotation_tasks = self.create_allele_balance_annotation_tasks(filter_variant_missingness_tasks, 5)
        # 6. annotate with 1000G
        annotate_1000G_tasks = self.create_1000G_annotation_tasks(allele_balance_annotation_tasks, 6)
        # 7. annotate with gnomAD
        annotate_gnomAD_tasks = self.create_gnomAD_annotation_tasks(annotate_1000G_tasks, 7)
        # 8. VEP annotation
        annotate_vep_tasks = self.create_vep_annotation_tasks(annotate_gnomAD_tasks, 8)
        # 9. CADD annotation
        annotate_cadd_tasks = self.create_cadd_annotation_tasks(annotate_vep_tasks, 9)
        # 10. Low-Confidence-Region annotation
        annotate_lcr_tasks = self.create_LCR_annotation_tasks(annotate_cadd_tasks, 10)
        # 11. LINSIGHT annotation
        annotate_linsight_tasks = self.create_LINSIGHT_annotation_tasks(annotate_lcr_tasks, 11)
        # 12. VCF concatenation
        concatenated_vcfs = self.create_concatenate_vcfs_task(annotate_linsight_tasks, 12)
        # 13. bcftools stats
        bcftools_stats_tasks = self.create_bcftools_stats_tasks(annotate_gnomAD_tasks, 13)
        # 13.1 Merge & Plot bcftools stats
        bcftools_stats_summary_task = self.create_bcftools_stats_summary_task(bcftools_stats_tasks, 13.1)
#        # 14. GATK VariantEval
#        variant_eval_tasks = self.create_variant_eval_tasks(annotate_gnomAD_tasks, 14)
#        # 14.1. Merge & Plot GATK VariantEval Stats
#        variant_eval_summary_task = self.create_variant_eval_summary_task(variant_eval_tasks, 14.1)

    def create_bcftools_stats_summary_task(self, parent_tasks, step_number):
        stage = self._construct_task_name('bcftools-stats-summary', step_number)
        output_dir = os.path.join(self.config.rootdir, stage)

        prior_stage_name = parent_tasks[0].stage.name
        input_dir = os.path.join(self.config.rootdir, prior_stage_name)

        lsf_params = get_lsf_params(
                bcftools_stats_summary_lsf_params,
                self.config.email,
                self.config.docker
        )
        lsf_params_json = to_json(lsf_params)

        task = {
            'func' : bcftools_stats_summary,
            'params' : {
                'in_dir' : input_dir,
                'out_dir' : output_dir,
            },
            'stage_name' : stage,
            'uid' : 'all-chroms',
            'drm_params' : lsf_params_json,
            'parents' : parent_tasks,
        }

        summary_task = self.workflow.add_task(**task)
        return summary_task

    def create_concatenate_vcfs_task(self, parent_tasks, step_number):
        tasks = list()
        stage = self._construct_task_name('concat-vcfs', step_number)
        output_dir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                concatenate_vcfs_lsf_params,
                self.config.email,
                self.config.docker,
                self.config.drm_queue
        )
        lsf_params_json = to_json(lsf_params)

        def region_key(task):
            reference_fai = '/gscmnt/gc2802/halllab/ccdg_resources/genomes/human/GRCh38DH/all_sequences.fa.fai'
            return Region(reference_fai, task.params['in_chrom'])

        def chromosome_key(task):
            reference_fai = '/gscmnt/gc2802/halllab/ccdg_resources/genomes/human/GRCh38DH/all_sequences.fa.fai'
            return Region(reference_fai, task.params['in_chrom']).chrom

        for ref_chrom, chrom_tasks in groupby(sorted(parent_tasks, key=region_key), key=chromosome_key):
            ptasks = list(chrom_tasks)
            input_vcfs = [ x.params['out_vcf'] for x in ptasks ]
            output_vcf = 'concatenated.c{}.vcf.gz'.format(ref_chrom)
            output_log = 'concatenate.{}.log'.format(ref_chrom)
            task = {
                'func' : concatenate_vcfs,
                'params' : {
                    'in_vcfs'  : input_vcfs,
                    'in_chrom' : ref_chrom,
                    'out_vcf' : os.path.join(output_dir, ref_chrom, output_vcf),
                    'out_log' : os.path.join(output_dir, ref_chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=ref_chrom),
                'drm_params' : lsf_params_json,
                'parents' : ptasks,
            }
            tasks.append( self.workflow.add_task(**task) )
        return tasks

    def create_variant_eval_summary_task(self, parent_tasks, step_number):
        stage = self._construct_task_name('gatk-variant-eval-summary', step_number)
        output_dir = os.path.join(self.config.rootdir, stage)

        prior_stage_name = parent_tasks[0].stage.name
        input_dir = os.path.join(self.config.rootdir, prior_stage_name)

        lsf_params = get_lsf_params(
                variant_eval_summary_lsf_params,
                self.config.email,
                self.config.docker
        )
        lsf_params_json = to_json(lsf_params)

        task = {
            'func' : variant_eval_summary,
            'params' : {
                'in_dir' : input_dir,
                'out_dir' : output_dir,
            },
            'stage_name' : stage,
            'uid' : 'all-chroms',
            'drm_params' : lsf_params_json,
            'parents' : parent_tasks,
        }

        summary_task = self.workflow.add_task(**task)
        return summary_task

    def create_bcftools_stats_tasks(self, parent_tasks, step_number):
        tasks = []
        stage = self._construct_task_name('bcftools-stats', step_number)
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                bcftools_stats_lsf_params,
                self.config.email,
                self.config.docker
        )
        lsf_params_json = to_json(lsf_params)

        for ptask in parent_tasks:
            chrom = ptask.params['in_chrom']
            output_stats = '{}.stats.out'.format(chrom)
            task = {
                'func' : bcftools_stats,
                'params' : {
                    'in_vcf' : ptask.params['out_vcf'],
                    'in_chrom' : chrom,
                    'out_stats' : os.path.join(basedir, chrom, output_stats),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' : lsf_params_json,
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_variant_eval_tasks(self, parent_tasks, step_number):
        tasks = []
        stage = self._construct_task_name('gatk-variant-eval', step_number)
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                gatk_variant_eval_lsf_params,
                self.config.email,
                self.config.docker
        )
        lsf_params_json = to_json(lsf_params)

        for ptask in parent_tasks:
            chrom = ptask.params['in_chrom']
            output_stats = 'chrom-{}-variant-eval.out'.format(chrom)
            output_log = 'chrom-{}-variant-eval.log'.format(chrom)
            task = {
                'func' : gatk_variant_eval,
                'params' : {
                    'in_vcf' : ptask.params['out_vcf'],
                    'in_chrom' : chrom,
                    'out_stats' : os.path.join(basedir, chrom, output_stats),
                    'out_log' : os.path.join(basedir, chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' : lsf_params_json,
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_LINSIGHT_annotation_tasks(self, parent_tasks, step_number):
        tasks = []
        stage = self._construct_task_name('LINSIGHT-annotation', step_number)
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                annotation_LINSIGHT_lsf_params,
                self.config.email,
                self.config.docker,
                self.config.drm_queue
        )
        lsf_params_json = to_json(lsf_params)

        for ptask in parent_tasks:
            chrom = ptask.params['in_chrom']
            output_vcf = 'b38.LINSIGHT.annotated.c{}.vcf.gz'.format(chrom)
            output_log = 'LINSIGHT.annotation.{}.log'.format(chrom)
            task = {
                'func' : annotation_LINSIGHT,
                'params' : {
                    'in_vcf'  : ptask.params['out_vcf'],
                    'in_chrom' : chrom,
                    'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                    'out_log' : os.path.join(basedir, chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' : lsf_params_json,
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_LCR_annotation_tasks(self, parent_tasks, step_number):
        tasks = []
        stage = self._construct_task_name('Low-Confidence-Region-annotation', step_number)
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                annotation_LCR_lsf_params,
                self.config.email,
                self.config.docker,
                self.config.drm_queue
        )
        lsf_params_json = to_json(lsf_params)

        for ptask in parent_tasks:
            chrom = ptask.params['in_chrom']
            output_vcf = 'b38.LCR.annotated.c{}.vcf.gz'.format(chrom)
            output_log = 'LCR.annotation.{}.log'.format(chrom)
            task = {
                'func' : annotation_LCR,
                'params' : {
                    'in_vcf'  : ptask.params['out_vcf'],
                    'in_chrom' : chrom,
                    'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                    'out_log' : os.path.join(basedir, chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' : lsf_params_json,
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_cadd_annotation_tasks(self, parent_tasks, step_number):
        tasks = []
        stage = self._construct_task_name('cadd-annotation', step_number)
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                annotation_cadd_lsf_params,
                self.config.email,
                self.config.docker,
                self.config.drm_queue
        )
        lsf_params_json = to_json(lsf_params)

        for ptask in parent_tasks:
            chrom = ptask.params['in_chrom']
            output_vcf = 'b38.cadd.annotated.c{}.vcf.gz'.format(chrom)
            output_log = 'cadd.annotation.{}.log'.format(chrom)
            task = {
                'func' : annotation_cadd,
                'params' : {
                    'in_vcf'  : ptask.params['out_vcf'],
                    'in_chrom' : chrom,
                    'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                    'out_log' : os.path.join(basedir, chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' : lsf_params_json,
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_vep_annotation_tasks(self, parent_tasks, step_number):
        tasks = []
        stage = self._construct_task_name('vep-annotation', step_number)
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                annotation_vep_lsf_params,
                self.config.email,
                self.config.docker,
                self.config.drm_queue
        )
        lsf_params_json = to_json(lsf_params)

        for ptask in parent_tasks:
            chrom = ptask.params['in_chrom']
            output_vcf = 'b38.vep.annotated.c{}.vcf.gz'.format(chrom)
            output_log = 'vep.annotation.{}.log'.format(chrom)
            task = {
                'func' : annotation_vep,
                'params' : {
                    'in_vcf'  : ptask.params['out_vcf'],
                    'in_chrom' : chrom,
                    'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                    'out_log' : os.path.join(basedir, chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' : lsf_params_json,
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_gnomAD_annotation_tasks(self, parent_tasks, step_number):
        tasks = []
        stage = self._construct_task_name('annotate-w-gnomAD', step_number)
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                annotation_gnomAD_lsf_params,
                self.config.email,
                self.config.docker
        )
        lsf_params_json = to_json(lsf_params)

        for ptask in parent_tasks:
            chrom = ptask.params['in_chrom']
            output_vcf = 'gnomAD-annotated.c{}.vcf.gz'.format(chrom)
            output_log = 'gnomAD-annotate.{}.log'.format(chrom)
            task = {
                'func' : annotation_gnomAD,
                'params' : {
                    'in_vcf' : ptask.params['out_vcf'],
                    'in_chrom' : chrom,
                    'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                    'out_log' : os.path.join(basedir, chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' : lsf_params_json,
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_1000G_annotation_tasks(self, parent_tasks, step_number):
        tasks = []
        stage = self._construct_task_name('annotate-w-1000G', step_number)
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                annotation_1000G_lsf_params,
                self.config.email,
                self.config.docker
        )
        lsf_params_json = to_json(lsf_params)

        for ptask in parent_tasks:
            chrom = ptask.params['in_chrom']
            output_vcf = '1kg-annotated.c{}.vcf.gz'.format(chrom)
            output_log = '1000G-annotate.{}.log'.format(chrom)
            task = {
                'func' : annotation_1000G,
                'params' : {
                    'in_vcf' : ptask.params['out_vcf'],
                    'in_chrom' : chrom,
                    'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                    'out_log' : os.path.join(basedir, chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' : lsf_params_json,
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_allele_balance_annotation_tasks(self, parent_tasks, step_number):
        tasks = []
        stage = self._construct_task_name('allele-balance-annotation', step_number)
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                annotate_allele_balances_lsf_params,
                self.config.email,
                self.config.docker,
                self.config.drm_queue
        )
        lsf_params_json = to_json(lsf_params)

        for ptask in parent_tasks:
            chrom = ptask.params['in_chrom']
            output_vcf = 'combined.c{chrom}.vcf.gz'.format(chrom=chrom)
            output_log = 'allele-balance-{}.log'.format(chrom)
            task = {
                'func' : annotate_allele_balances,
                'params' : {
                    'in_vcf' : ptask.params['out_vcf'],
                    'in_chrom' : chrom,
                    'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                    'out_log' : os.path.join(basedir, chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' : lsf_params_json,
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_filter_variant_missingness_tasks(self, parent_tasks, step_number):
        tasks = []
        stage = self._construct_task_name('filter-variant-missingness', step_number)
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                filter_variant_missingness_lsf_params,
                self.config.email,
                self.config.docker
        )
        lsf_params_json = to_json(lsf_params)

        for ptask in parent_tasks:
            chrom = ptask.params['in_chrom']
            output_vcf = 'combined.c{chrom}.vcf.gz'.format(chrom=chrom)
            output_log = 'filter-missingness-{}.log'.format(chrom)
            task = {
                'func' : filter_variant_missingness,
                'params' : {
                    'in_vcf' : ptask.params['out_vcf'],
                    'in_chrom' : chrom,
                    'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                    'out_log' : os.path.join(basedir, chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' : lsf_params_json,
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_remove_symbolic_deletion_tasks(self, parent_tasks, step_number):
        tasks = []
        stage = self._construct_task_name('remove-symbolic-alleles', step_number)
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                remove_symbolic_deletion_alleles_lsf_params,
                self.config.email,
                self.config.docker
                )
        lsf_params_json = to_json(lsf_params)

        for ptask in parent_tasks:
            chrom = ptask.params['in_chrom']
            output_vcf = 'combined.c{chrom}.vcf.gz'.format(chrom=chrom)
            output_log = 'remove-symbolic-alleles-chrom-{}.log'.format(chrom)
            task = {
                    'func'   : remove_symbolic_deletion_alleles,
                    'params' : {
                        'in_vcf' : ptask.params['out_vcf'],
                        'in_chrom' : chrom,
                        'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                        'out_log' : os.path.join(basedir, chrom, output_log),
                        },
                    'stage_name' : stage,
                    'uid' : '{chrom}'.format(chrom=chrom),
                    'drm_params' : lsf_params_json,
                    'parents' : [ptask],
                    }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_decompose_normalize_unique_tasks(self, step_number):
        tasks = []
        stage = self._construct_task_name('decompose-normalize-uniq', step_number)
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                normalize_decompose_unique_lsf_params,
                self.config.email,
                self.config.docker
        )
        lsf_params_json = to_json(lsf_params)

        for chrom in self.config.chroms:
            output_vcf = 'combined.c{chrom}.vcf.gz'.format(chrom=chrom)
            output_log = 'decompose-normalize-unique-{}.log'.format(chrom)
            task = {
                'func' : normalize_decompose_unique,
                'params' : {
                    'in_vcf' : self.config.vcfs[chrom],
                    'in_chrom' : chrom,
                    'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                    'out_log' : os.path.join(basedir, chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' : lsf_params_json,
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_calculate_sample_missingness_task(self, parent_tasks, step_number):
        stage = self._construct_task_name('calculate-sample-missingness', step_number)
        output_dir = os.path.join(self.config.rootdir, stage)

        prior_stage_name = parent_tasks[0].stage.name
        input_dir = os.path.join(self.config.rootdir, prior_stage_name)
        input_json_wildcard_path = os.path.join(input_dir, '*', '*.json')

        lsf_params = get_lsf_params(
                calculate_sample_missingness_lsf_params,
                self.config.email,
                self.config.docker
        )
        lsf_params_json = to_json(lsf_params)

        task = {
            'func' : calculate_sample_missingness,
            'params' : {
                'in_json' : input_json_wildcard_path,
                'out_stats' : os.path.join(output_dir, 'sample-missingness-pct.dat'),
                'out_log' : os.path.join(output_dir, 'sample-missingness-pct.dat.log'),
            },
            'stage_name' : stage,
            'uid' : '1-22',
            'drm_params' : lsf_params_json,
            'parents' : parent_tasks,
        }

        summary_task = self.workflow.add_task(**task)
        return summary_task

    def create_count_sample_missingness_tasks(self, step_number):
        tasks = []
        stage = self._construct_task_name('count-sample-missingness', step_number)
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                count_sample_missingness_lsf_params,
                self.config.email,
                self.config.docker
        )
        lsf_params_json = to_json(lsf_params)

        for chrom in self.config.chroms:

            # only count missing genotypes on chromosomes 1-22 (not X, Y, or MT)
            if not chrom[0].isdigit() : continue

            output_json = '{chrom}-sample-missingness-counts.json'.format(chrom=chrom)
            output_log = '{}-sample-missingness-counts.log'.format(chrom)
            task = {
                'func' : count_sample_missingness,
                'params' : {
                    'in_vcf' : self.config.vcfs[chrom],
                    'in_chrom' : chrom,
                    'out_json' : os.path.join(basedir, chrom, output_json),
                    'out_log' : os.path.join(basedir, chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' : lsf_params_json,
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def _construct_task_name(self, name, number):
        task_name = '{}-{}'.format(number, name)
        return task_name

# C M D S #####################################################################
def get_lsf_params(task_lsf_fn, email, docker, queue):
    lsf_params = task_lsf_fn(email, queue)
    if docker and (lsf_params['q'] != 'research-hpc'):
        lsf_params['a'] = "'docker(registry.gsc.wustl.edu/genome/genome_perl_environment:23)'"
        lsf_params['q'] = "research-hpc"
        lsf_params['M'] = 16000000
        lsf_params['R'] = 'select[mem>10000 && ncpus>8] rusage[mem=16000]'
    return lsf_params

def bcftools_stats_summary(in_dir, out_dir):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/bcftools-stats-summary-plots.sh'),
    }
    cmd_args = merge_params(default, args)
    cmd = "{script} {in_dir} {out_dir}".format(**cmd_args)
    return cmd

def bcftools_stats_summary_lsf_params(email, queue):
    return  {
        'u' : email,
        'N' : None,
        'q' : queue,
        'M' : 8000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }

def bcftools_stats(in_vcf, in_chrom, out_stats):
    args = locals()
    default = {
        'bcftools' : '/gscmnt/gc2802/halllab/idas/software/local/bin/bcftools1.4',
        'reference' : '/gscmnt/gc2802/halllab/ccdg_resources/genomes/human/GRCh38DH/all_sequences.fa',
    }

    cmd_args = merge_params(default, args)

    cmd = ( "{bcftools} stats "
            "--split-by-ID "
            "-F {reference} "
            "-s - "
            "-f '.,PASS' "
            "{in_vcf} "
            ">{out_stats}").format(**cmd_args)

    return cmd

def bcftools_stats_lsf_params(email, queue):
    return  {
        'u' : email,
        'N' : None,
        'q' : queue,
        'M' : 10000000,
        'R' : 'select[mem>10000 && ncpus>8] rusage[mem=10000]',
    }

def concatenate_vcfs(in_vcfs, in_chrom, out_vcf, out_log):
    args = locals()
    default = {
        'bcftools' : '/gscmnt/gc2802/halllab/idas/software/local/bin/bcftools1.4',
    }

    cmd_args = merge_params(default, args)
    cmd = None
    if len(cmd_args['in_vcfs']) == 1:
        cmd_args['in_vcfs'] = cmd_args['in_vcfs'][0]
        cmd = ( "cp -v {in_vcfs} {out_vcf}"
                ">{out_log} "
                "2>&1").format(**cmd_args)
    else:
        cmd_args['in_vcfs'] = ' '.join([ x for x in in_vcfs if not empty_gzipped_vcf(x) ])
        cmd = ( "{bcftools} concat "
                "-a "
                "{in_vcfs} "
                "-O z -o {out_vcf}"
                ">{out_log} "
                "2>&1").format(**cmd_args)
    return cmd

def concatenate_vcfs_lsf_params(email, queue):
    return  {
        'u' : email,
        'N' : None,
        'q' : queue,
        'M' : 10000000,
        'R' : 'select[mem>10000 && ncpus>8] rusage[mem=10000]',
    }

def variant_eval_summary(in_dir, out_dir):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/merge-and-plot-gatk-variant-eval-stats.sh'),
    }
    cmd_args = merge_params(default, args)
    cmd = "{script} {in_dir} {out_dir}".format(**cmd_args)
    return cmd

def variant_eval_summary_lsf_params(email, queue):
    return  {
        'u' : email,
        'N' : None,
        'q' : queue,
        'M' : 8000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }

def gatk_variant_eval(in_chrom, in_vcf, out_stats, out_log):
    args = locals()
    default = {
        'java' : '/gapp/x64linux/opt/java/jre/jre1.8.0_31/bin/java',
        'jar'  : '/gscmnt/gc2802/halllab/idas/jira/BIO-1662/vendor/local/jars/GenomeAnalysisTK-3.5-idas-experimental-293f64d-2016.02.19.jar',
        'java_opts' : "-Xmx4096m",
        'reference' : '/gscmnt/gc2719/halllab/genomes/human/GRCh37/1kg_phase1/human_g1k_v37.fasta',
        'dbsnp': '/gscmnt/gc2802/halllab/idas/jira/BIO-1662/data/derived/FinnMetSeq-WGS/10-decompose-normalize-1000G-variant-ref-v1/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.decompose.normalize.vcf.gz',
    }

    cmd_args = merge_params(default, args)

    cmd = ( "{java} -jar {jar} "
            "-nt 8 "
            "-T VariantEval "
            "-D {dbsnp} "
            "-R {reference} "
            "-ST Sample "
            "-noST "
            "-EV TiTvVariantEvaluator "
            "-EV CountVariants "
            "-EV CompOverlap "
            "-EV IndelSummary "
            "-EV MultiallelicSummary "
            "-noEV "
            "-L {in_chrom} "
            "-eval {in_vcf} "
            "-o {out_stats} "
            ">{out_log} "
            "2>&1").format(**cmd_args)

    return cmd

def gatk_variant_eval_lsf_params(email, queue):
    return  {
        'u' : email,
        'N' : None,
        'q' : queue,
        'M' : 10000000,
        'R' : 'select[mem>10000 && ncpus>8] rusage[mem=10000]',
    }

def annotation_LINSIGHT(in_vcf, in_chrom, out_vcf, out_log):
    args = locals()
    default = {
        'main_script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/annotate-w-LINSIGHT.sh'),
    }
    cmd_args = merge_params(default, args)
    cmd = ("{main_script} "
           "{in_vcf} "
           "{out_vcf} "
           ">{out_log} 2>&1" ).format(**cmd_args)
    return cmd

def annotation_LINSIGHT_lsf_params(email, queue):
    return  {
        'u' : email,
        'N' : None,
        'q' : queue,
        'M' : 8000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }

def annotation_LCR(in_vcf, in_chrom, out_vcf, out_log):
    args = locals()
    default = {
        'main_script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/annotate-regions-of-low-confidence.sh'),
    }
    cmd_args = merge_params(default, args)
    cmd = ("{main_script} "
           "{in_vcf} "
           "{out_vcf} "
           ">{out_log} 2>&1" ).format(**cmd_args)
    return cmd

def annotation_LCR_lsf_params(email, queue):
    return  {
        'u' : email,
        'N' : None,
        'q' : queue,
        'M' : 8000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }

def annotation_cadd(in_vcf, in_chrom, out_vcf, out_log):
    args = locals()
    default = {
        'main_script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/run-cadd.sh'),
        'merge_script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/merge-in-cadd.py'),
        'b37_to_b38_integration_script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/integrate-b37-annotations-to-b38.py'),
    }
    cmd_args = merge_params(default, args)
    cmd = ("{main_script} "
           "{in_vcf} "
           "{out_vcf} "
           "{merge_script} "
           "{b37_to_b38_integration_script} "
           ">{out_log} 2>&1" ).format(**cmd_args)
    return cmd

def annotation_cadd_lsf_params(email, queue):
    return  {
        'u' : email,
        'N' : None,
        'q' : queue,
        'M' : 64000000,
        'R' : 'select[mem>60000 && ncpus>8] rusage[mem=64000]',
    }

def annotation_vep(in_vcf, in_chrom, out_vcf, out_log):
    args = locals()
    default = {
        'main_script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/run-vep.sh'),
    }
    cmd_args = merge_params(default, args)
    cmd = ("{main_script} "
           "{in_vcf} "
           "{out_vcf} "
           ">{out_log} 2>&1" ).format(**cmd_args)
    return cmd

def annotation_vep_lsf_params(email, queue):
    return  {
        'u' : email,
        'N' : None,
        'a' : "'docker(willmclaren/ensembl-vep:release_88)'",
        'q' : 'research-hpc',
        'M' : 68000000,
        'R' : 'select[mem>60000 && ncpus>8] rusage[mem=68000]',
    }

def annotation_gnomAD(in_vcf, in_chrom, out_vcf, out_log):
    args = locals()
    default = {
        'main_script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/annotate-w-gnomAD.sh'),
        'b37_to_b38_integration_script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/integrate-b37-annotations-to-b38.py'),
    }
    cmd_args = merge_params(default, args)
    cmd = ("{main_script} "
           "{in_chrom} "
           "{in_vcf} "
           "{out_vcf} "
           "{b37_to_b38_integration_script} "
           ">{out_log} 2>&1" ).format(**cmd_args)
    return cmd

def annotation_gnomAD_lsf_params(email, queue):
    return  {
        'u' : email,
        'N' : None,
        'q' : queue,
        'M' : 8000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }

def annotation_1000G(in_vcf, in_chrom, out_vcf, out_log):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/annotate-w-1000G.sh'),
        'integrate_script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/integrate-b37-annotations-to-b38.py'),
    }
    cmd_args = merge_params(default, args)
    cmd = "{script} {in_vcf} {out_vcf} {integrate_script} >{out_log} 2>&1".format(**cmd_args)
    return cmd

def annotation_1000G_lsf_params(email, queue):
    return  {
        'u' : email,
        'N' : None,
        'q' : queue,
        'M' : 8000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }

def annotate_allele_balances(in_vcf, in_chrom, out_vcf, out_log):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr/allele-balance-annotation.sh'),
        'python_script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr/annotate-allele-balances.py'),
        'python_executable' : sys.executable,
    }
    cmd_args = merge_params(default, args)

    cmd = ( "{script} "
            "{python_executable} {python_script} "
            "{in_vcf} {out_vcf} {in_chrom} "
            ">{out_log} 2>&1" ).format(**cmd_args)

    return cmd

def annotate_allele_balances_lsf_params(email, queue):
    return  {
        'u' : email,
        'N' : None,
        'q' : queue,
        'M' : 8000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }

def filter_variant_missingness(in_vcf, in_chrom, out_vcf, out_log):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/filter-missingness-sites.sh'),
        'python_script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/filter-site-missingness.py'),
        'python_executable' : sys.executable,
    }
    cmd_args = merge_params(default, args)

    if in_chrom.startswith('Y') or in_chrom.startswith('y'):
        cmd_args['out_vcf'] = os.path.dirname(out_vcf)
        cmd = "/bin/cp -v {in_vcf}* {out_vcf} >{out_log} 2>&1".format(**cmd_args)
    else:
        cmd = ( "{script} "
                "{python_executable} {python_script} "
                "{in_vcf} {out_vcf} {in_chrom} "
                ">{out_log} 2>&1" ).format(**cmd_args)
    return cmd

def filter_variant_missingness_lsf_params(email, queue):
    return  {
        'u' : email,
        'N' : None,
        'q' : queue,
        'M' : 8000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }

def remove_symbolic_deletion_alleles(in_vcf, in_chrom, out_vcf, out_log):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/remove-symbolic.sh'),
    }

    cmd_args = merge_params(default, args)

    cmd = "{script} {in_vcf} {out_vcf} >{out_log} 2>&1".format(**cmd_args)
    return cmd

def remove_symbolic_deletion_alleles_lsf_params(email, queue):
    return {
        'u' : email,
        'N' : None,
        'q' : queue,
        'M' : 4000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }

def normalize_decompose_unique(in_vcf, in_chrom, out_vcf, out_log):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/run-decompose.sh'),
    }
    cmd_args = merge_params(default, args)
    cmd = "{script} {in_vcf} {out_vcf} {in_chrom} >{out_log} 2>&1".format(**cmd_args)
    return cmd

def normalize_decompose_unique_lsf_params(email, queue):
    return  {
        'u' : email,
        'N' : None,
        'q' : queue,
        'M' : 32000000,
        'R' : 'select[mem>32000 && ncpus>8] rusage[mem=32000]',
    }

def calculate_sample_missingness(in_json, out_stats, out_log):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/calculate-overall-sample-missingness.py'),
        'python' : sys.executable,
    }
    cmd_args = merge_params(default, args)
    cmd = "{python} {script} --out={out_stats} {in_json} >{out_log} 2>&1".format(**cmd_args)
    return cmd

def calculate_sample_missingness_lsf_params(email, queue):
    return  {
        'u' : email,
        'N' : None,
        'q' : queue,
        'M' : 8000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }

def count_sample_missingness(in_vcf, in_chrom, out_json, out_log):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/count-sample-missingness.py'),
        'python' : sys.executable,
    }
    cmd_args = merge_params(default, args)
    cmd = "{python} {script} --out={out_json} {in_vcf} >{out_log} 2>&1".format(**cmd_args)
    return cmd

def count_sample_missingness_lsf_params(email, queue):
    return  {
        'u' : email,
        'N' : None,
        'q' : queue,
        'M' : 8000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }
