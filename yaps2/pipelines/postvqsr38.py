import os, pwd, sys
import pkg_resources
from cosmos.api import Cosmos, Dependency, default_get_submit_args
from yaps2.utils import to_json, merge_params, natural_key

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
        # 1. remove unused alternates
        remove_ac_0_tasks = self.create_remove_ac_0_tasks()
        # 2. calculate sample missingness (counting phase)
        count_sample_missingness_tasks = self.create_count_sample_missingness_tasks(remove_ac_0_tasks)
        # 2.1 calculate sample missingness (merge and calculation phase)
        calculate_sample_missingness_task = self.create_calculate_sample_missingness_task(count_sample_missingness_tasks)
        # 3. denormalize, decompose, and uniq
        dnu_tasks = self.create_decompose_normalize_unique_tasks(remove_ac_0_tasks)
        # 4. remove symbolic alleles
        rsa_tasks = self.create_remove_symbolic_deletion_tasks(dnu_tasks)
        # 5. filter missingess
        filter_variant_missingness_tasks = self.create_filter_variant_missingness_tasks(rsa_tasks)
        # 6. annotate with 1000G
        annotate_1000G_tasks = self.create_1000G_annotation_tasks(filter_variant_missingness_tasks)
        # 7. annotate with ExAC
        annotate_ExAC_tasks = self.create_ExAC_annotation_tasks(annotate_1000G_tasks)
        # 8. CADD/VEP annotation
        # annotate_vep_cadd_task = self.create_vep_cadd_annotation_task(annotate_ExAC_tasks)
        # 9. GATK VariantEval
        variant_eval_tasks = self.create_variant_eval_tasks(annotate_ExAC_tasks)
        # 9.1. Merge & Plot GATK VariantEval Stats
        variant_eval_summary_task = self.create_variant_eval_summary_task(variant_eval_tasks)
        # 10. bcftools stats
        bcftools_stats_tasks = self.create_bcftools_stats_tasks(annotate_ExAC_tasks)
        # 10.1 Merge & Plot bcftools stats
        bcftools_stats_summary_task = self.create_bcftools_stats_summary_task(bcftools_stats_tasks)

    def create_bcftools_stats_summary_task(self, parent_tasks):
        stage = '10.1-bcftools-stats-summary'
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

    def create_variant_eval_summary_task(self, parent_tasks):
        stage = '9.1-gatk-variant-eval-summary'
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

    def create_bcftools_stats_tasks(self, parent_tasks):
        tasks = []
        stage = '10-bcftools-stats'
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

    def create_variant_eval_tasks(self, parent_tasks):
        tasks = []
        stage = '9-gatk-variant-eval'
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

    def create_vep_cadd_annotation_task(self, parent_tasks):
        stage = '8-vep-cadd-annotation'
        output_dir = os.path.join(self.config.rootdir, stage)

        prior_stage_name = parent_tasks[0].stage.name
        input_dir = os.path.join(self.config.rootdir, prior_stage_name)

        lsf_params = get_lsf_params(
                annotation_VEP_CADD_lsf_params,
                self.config.email,
                self.config.docker
        )
        lsf_params_json = to_json(lsf_params)

        task = {
            'func' : annotation_VEP_CADD,
            'params' : {
                'in_dir' : input_dir,
                'out_dir' : output_dir,
            },
            'stage_name' : stage,
            'uid' : 'all-chroms',
            'drm_params' : lsf_params_json,
            'parents' : parent_tasks,
        }

        vep_cadd_task = self.workflow.add_task(**task)
        return vep_cadd_task

    def create_ExAC_annotation_tasks(self, parent_tasks):
        tasks = []
        stage = '7-annotate-w-ExAC'
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                annotation_ExAC_lsf_params,
                self.config.email,
                self.config.docker
        )
        lsf_params_json = to_json(lsf_params)

        for ptask in parent_tasks:
            chrom = ptask.params['in_chrom']
            output_vcf = 'ExAC-annotated.c{}.vcf.gz'.format(chrom)
            output_log = 'ExAC-annotate.{}.log'.format(chrom)
            task = {
                'func' : annotation_ExAC,
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

    def create_1000G_annotation_tasks(self, parent_tasks):
        tasks = []
        stage = '6-annotate-w-1000G'
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

    def create_filter_variant_missingness_tasks(self, parent_tasks):
        tasks = []
        stage = '5-filter-variant-missingness'
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
            output_stats = '{chrom}.stats.missingness.out'.format(chrom=chrom)
            output_log = 'filter-missingness-{}.log'.format(chrom)
            task = {
                'func' : filter_variant_missingness,
                'params' : {
                    'in_vcf' : ptask.params['out_vcf'],
                    'in_chrom' : chrom,
                    'out_vcf' : os.path.join(basedir, chrom, output_vcf),
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

    def create_decompose_normalize_unique_tasks(self, parent_tasks):
        tasks = []
        stage = '3-decompose-normalize-uniq'
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                normalize_decompose_unique_lsf_params,
                self.config.email,
                self.config.docker
        )
        lsf_params_json = to_json(lsf_params)

        for ptask in parent_tasks:
            chrom = ptask.params['in_chrom']
            output_vcf = 'combined.c{chrom}.vcf.gz'.format(chrom=chrom)
            output_log = 'decompose-normalize-unique-{}.log'.format(chrom)
            task = {
                'func' : normalize_decompose_unique,
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

    def create_calculate_sample_missingness_task(self, parent_tasks):
        stage = '2.1-calculate-sample-missingness'
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

    def create_count_sample_missingness_tasks(self, parent_tasks):
        tasks = []
        stage = '2-count-sample-missingness'
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                count_sample_missingness_lsf_params,
                self.config.email,
                self.config.docker
        )
        lsf_params_json = to_json(lsf_params)

        for ptask in parent_tasks:
            chrom = ptask.params['in_chrom']

            # only count missing genotypes on chromosomes 1-22 (not X, Y, or MT)
            if not chrom[0].isdigit() : continue

            output_json = '{chrom}-sample-missingness-counts.json'.format(chrom=chrom)
            output_log = '{}-sample-missingness-counts.log'.format(chrom)
            task = {
                'func' : count_sample_missingness,
                'params' : {
                    'in_vcf' : ptask.params['out_vcf'],
                    'in_chrom' : chrom,
                    'out_json' : os.path.join(basedir, chrom, output_json),
                    'out_log' : os.path.join(basedir, chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' : lsf_params_json,
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_remove_ac_0_tasks(self):
        tasks = []
        stage = '1-select-variants-ac-0-removal'
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                gatk_select_variants_remove_ac_0_lsf_params,
                self.config.email,
                self.config.docker
        )
        lsf_params_json = to_json(lsf_params)

        for chrom in self.config.chroms:
            vcf = self.config.vcfs[chrom]
            output_vcf = 'combined.c{chrom}.vcf.gz'.format(chrom=chrom)
            output_log = 'select-variants-chrom-{}-gatk.log'.format(chrom)
            task = {
                'func'   : gatk_select_variants_remove_ac_0,
                'params' : {
                    'in_chrom' : chrom,
                    'in_vcf' : vcf,
                    'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                    'out_log' : os.path.join(basedir, chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' : lsf_params_json,
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_remove_symbolic_deletion_tasks(self):
        tasks = []
        stage = '4-remove-symbolic-alleles'
        basedir = os.path.join(self.config.rootdir, stage)

        lsf_params = get_lsf_params(
                remove_symbolic_deletion_alleles_lsf_params,
                self.config.email,
                self.config.docker
                )
        lsf_params_json = to_json(lsf_params)

        for chrom in self.config.chroms:
            vcf = self.config.vcfs[chrom]
            output_vcf = 'combined.c{chrom}.vcf.gz'.format(chrom=chrom)
            output_log = 'remove-symbolic-alleles-chrom-{}-gatk.log'.format(chrom)
            task = {
                    'func'   : remove_symbolic_deletion_alleles,
                    'params' : {
                        'in_vcf' : vcf,
                        'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                        'out_log' : os.path.join(basedir, chrom, output_log),
                        },
                    'stage_name' : stage,
                    'uid' : '{chrom}'.format(chrom=chrom),
                    'drm_params' : lsf_params_json,
                    }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

# C M D S #####################################################################
def get_lsf_params(task_lsf_fn, email, docker):
    lsf_params = task_lsf_fn(email)
    if docker:
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

def bcftools_stats_summary_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "ccdg",
        'M' : 8000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }

def bcftools_stats(in_vcf, in_chrom, out_stats):
    args = locals()
    default = {
        'bcftools' : '/gsc/bin/bcftools1.2',
        'reference' : '/gscmnt/ams1102/info/model_data/2869585698/build106942997/all_sequences.fa',
    }

    cmd_args = merge_params(default, args)

    cmd = ( "{bcftools} stats "
            "--split-by-ID "
            "-F {reference} "
            "-f '.,PASS' "
            "{in_vcf} "
            ">{out_stats}").format(**cmd_args)

    return cmd

def bcftools_stats_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "ccdg",
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

def variant_eval_summary_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "ccdg",
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

def gatk_variant_eval_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "ccdg",
        'M' : 10000000,
        'R' : 'select[mem>10000 && ncpus>8] rusage[mem=10000]',
    }

def annotation_VEP_CADD(in_dir, out_dir):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/vep-cadd-annotation.sh'),
    }
    cmd_args = merge_params(default, args)
    cmd = "{script} {in_dir} {out_dir}".format(**cmd_args)
    return cmd

def annotation_VEP_CADD_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "ccdg",
        'M' : 8000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }

def annotation_ExAC(in_vcf, in_chrom, out_vcf, out_log):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/annotate-w-ExAC.sh'),
    }
    cmd_args = merge_params(default, args)
    cmd = "{script} {in_vcf} {out_vcf} >{out_log} 2>&1".format(**cmd_args)
    return cmd

def annotation_ExAC_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "ccdg",
        'M' : 8000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }

def annotation_1000G(in_vcf, in_chrom, out_vcf, out_log):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/annotate-w-1000G.sh'),
    }
    cmd_args = merge_params(default, args)
    cmd = "{script} {in_vcf} {out_vcf} >{out_log} 2>&1".format(**cmd_args)
    return cmd

def annotation_1000G_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "ccdg",
        'M' : 8000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }

def filter_variant_missingness(in_vcf, in_chrom, out_vcf, out_stats, out_log):
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
                "{in_vcf} {out_vcf} {out_stats} {in_chrom} "
                ">{out_log} 2>&1" ).format(**cmd_args)
    return cmd

def filter_variant_missingness_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "ccdg",
        'M' : 8000000,
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

def normalize_decompose_unique_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "ccdg",
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

def calculate_sample_missingness_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "ccdg",
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

def count_sample_missingness_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "ccdg",
        'M' : 8000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }

def gatk_select_variants_remove_ac_0(in_chrom, in_vcf, out_vcf, out_log):
    args = locals()
    default = {
        'java' : '/gapp/x64linux/opt/java/jre/jre1.7.0_45/bin/java',
        'jar'  : '/usr/share/java/GenomeAnalysisTK-3.4.jar',
        'java_opts' : "-Xmx4096m",
        'reference' : '/gscmnt/gc2719/halllab/genomes/human/GRCh37/1kg_phase1/human_g1k_v37.fasta',
    }

    cmd_args = merge_params(default, args)

    cmd = ( "{java} -jar {jar} "
            "-T SelectVariants -R {reference} "
            "--removeUnusedAlternates "
            "-V {in_vcf} "
            "-L {in_chrom} "
            "-o {out_vcf} "
            ">{out_log} "
            "2>&1").format(**cmd_args)

    return cmd

def gatk_select_variants_remove_ac_0_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "ccdg",
        'M' : 8000000,
        'R' : 'select[mem>8000 && ncpus>8] rusage[mem=8000]',
    }

def remove_symbolic_deletion_alleles(in_vcf, out_vcf, out_log):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr38/remove-symbolic.sh'),
    }

    cmd_args = merge_params(default, args)

    cmd = "{script} {in_vcf} {out_vcf} >{out_log} 2>&1".format(**cmd_args)
    return cmd

def remove_symbolic_deletion_alleles_lsf_params(email):
    return {
        'u' : email,
        'N' : None,
        'q' : "ccdg",
        'M' : 4000000,
        'R' : 'select[mem>8000] && ncpus>8] rusage[mem=8000]',
    }

