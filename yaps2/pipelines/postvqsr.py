import os, pwd
import pkg_resources
from cosmos.api import Cosmos, Dependency, default_get_submit_args
from yaps2.utils import to_json, merge_params, natural_key

class Config(object):
    def __init__(job_db, input_vcf_list, project_name, email, workspace):
        self.email = email
        self.db = job_db
        self.project_name = project_name
        self.rootdir = workspace

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
    def __init__(config, drm, restart):
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
        )

        self.setup_pipeline()

    def setup_pipeline(self):
        self.construct_pipeline()
        self.workflow.make_output_dirs()

    def run(self):
	# put set_successful to False if you intend to add more tasks to the
	# pipeline later
        self.workflow.run(set_successful=False)

    def construct_pipeline(self):
        # 1. remove unused alternates
        # 2. denormalize, decompose, and uniq
        remove_ac_0_tasks = self.create_remove_ac_0_tasks()
        dnu_tasks = self.create_decompose_normalize_unique_tasks(remove_ac_0_tasks)
        filter_missingness_tasks = self.create_filter_missingness_tasks(dnu_tasks)
        annotate_1000G_tasks = self.create_1000G_annotation_tasks(filter_missingness_tasks)

    def create_1000G_annotation_tasks(self, parent_tasks):
        tasks = []
        stage = '4-annotate-w-1000G'
        basedir = os.path.join(self.config.rootdir, stage)
        email = self.config.email

        for ptask in parent_tasks:
            chrom = ptask.params['in_chrom']
            output_vcf = '1kg-annotated.c{}.vcf.gz'.format(chrom)
            output_log = '1000G-annotate.{}.log'.format(chrom)
            task = {
                'func' : filter_missingness,
                'params' : {
                    'in_vcf' : ptask.params['out_vcf'],
                    'in_chrom' : chrom,
                    'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                    'out_log' : os.path.join(basedir, chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' :
                    to_json(annotation_1000G_lsf_params(email)),
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_filter_missingness_tasks(self, parent_tasks):
        tasks = []
        stage = '3-filter-missingness'
        basedir = os.path.join(self.config.rootdir, stage)
        email = self.config.email

        for ptask in parent_tasks:
            chrom = ptask.params['in_chrom']
            output_vcf = 'combined.c{chrom}.vcf.gz'.format(chrom=chrom)
            output_log = 'filter-missingness-{}.log'.format(chrom)
            task = {
                'func' : filter_missingness,
                'params' : {
                    'in_vcf' : ptask.params['out_vcf'],
                    'in_chrom' : chrom,
                    'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                    'out_log' : os.path.join(basedir, chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' :
                    to_json(normalize_decompose_unique_lsf_params(email)),
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_decompose_normalize_unique_tasks(self, parent_tasks):
        tasks = []
        stage = '2-decompose-normalize-uniq'
        basedir = os.path.join(self.config.rootdir, stage)
        email = self.config.email

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
                'drm_params' :
                    to_json(normalize_decompose_unique_lsf_params(email)),
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_remove_ac_0_tasks(self):
        tasks = []
        stage = '1-select-variants-ac-0-removal'
        basedir = os.path.join(self.config.rootdir, stage)
        email = self.config.email

        for chrom in self.config.chroms:
            vcf = self.config.vcfs[chrom]
            output_vcf = 'combined.c{chrom}.vcf.gz'.format(chrom=chrom)
            output_log = 'select-variants-chrom-{}-gatk.log'.format(chrom)
            task = {
                'func'   : gatk_select_variants_remove_ac_0,
                'params' : {
                    'in_chrom' : chrom
                    'in_vcf' : vcf,
                    'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                    'out_log' : os.path.join(basedir, chrom, output_log),
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' :
                    to_json(gatk_select_variants_remove_ac_0_lsf_params(email)),
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

# C M D S #####################################################################
def annotation_1000G(in_vcf, in_chrom, out_vcf, out_log):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr/annotate-with-1000G.sh'),
    }
    cmd_args = merge_params(default, args)
    cmd = "{script} {in_vcf} {out_vcf} 2>&1 >{out_log}".format(**cmd_args)
    return cmd

def annotation_1000G_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "long",
        'M' : 8000000,
        'R' : 'select[mem>8000] rusage[mem=8000]',
    }

def normalize_decompose_unique(in_vcf, in_chrom, out_vcf, out_log):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/postvqsr/run-decompose.sh'),
    }
    cmd_args = merge_params(default, args)
    cmd = "{script} {in_vcf} {out_vcf} 2>&1 >{out_log}".format(**cmd_args)
    return cmd

def normalize_decompose_unique_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "long",
        'M' : 8000000,
        'R' : 'select[mem>8000] rusage[mem=8000]',
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
            "-V {invcf} "
            "-L {chrom} "
            "-o {outvcf} "
            "2>&1 "
            ">{out_log}").format(**cmd_args)

    return cmd

def gatk_select_variants_remove_ac_0_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "long",
        'M' : 8000000,
        'R' : 'select[mem>8000] rusage[mem=8000]',
    }
