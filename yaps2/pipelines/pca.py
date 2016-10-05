from __future__ import print_function, division

import os, sys, pwd
import pkg_resources
from cosmos.api import Cosmos, Dependency, default_get_submit_args
from yaps2.utils import to_json, merge_params, natural_key, ensure_directory

class Config(object):
    def __init__(self, job_db, 
                 input_vcfs_file, project_name, email, workspace, vqslod_threshold):
        self.email = email
        self.db = job_db
        self.project_name = project_name
        self.rootdir = workspace
        self.vqslod_threshold = vqslod_threshold

        self.ensure_rootdir()

        if self.email is None:
            self.email = self.setup_email()

        if self.db is None:
            self.db = os.path.join(
                os.path.abspath(self.rootdir),
                '.job_queue.db'
            )

        self.vcfs = self.collect_input_vcfs(input_vcfs_file)
        self.chroms = self.get_ordered_chroms()

        # ensure that the sex chroms are excluded from analysis
        sex_chroms = ('x', 'y')
        for chrom in sex_chroms:
            self.ensure_chromosome_not_present(chrom, input_vcfs_file)

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

    def ensure_chromosome_not_present(self, chrom, input_file):
        lower_chroms = [ i.lower() for i in self.chroms ]
        if chrom in lower_chroms:
            msg = "Please exclude chromosome 'X' from {}".format(input_file)
            sys.exit(msg)

class Pipeline(object):
    def __init__(self, config, drm, restart):
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
        custom_log_dir = lambda task : os.path.join(self.config.rootdir, 'logs', task.stage.name, task.uid)
        self.workflow.run(set_successful=False, log_out_dir_func=custom_log_dir)

    def construct_pipeline(self):
        filter_biallelic_snps_tasks = self.create_filter_biallelic_snps_tasks()
        plink_binary_tasks = self.create_plink_binary_tasks(filter_biallelic_snps_tasks)
        plink_ld_prune_tasks = self.create_plink_ld_prune_tasks(plink_binary_tasks)
        plink_extract_prune_tasks = self.create_plink_extract_prune_tasks(plink_ld_prune_tasks)
        plink_merge_prune_files_task = self.create_plink_merge_prune_file_task(plink_extract_prune_tasks)
        eigenstrat_task = self.create_eigenstrat_smartpca_task(plink_merge_prune_files_task)
        data_frame_task = self.create_data_frame_task(eigenstrat_task)

    def create_data_frame_task(self, parent_task):
        stage = '7-make-data-frame'
        basedir = os.path.join(self.config.rootdir, stage)
        email = self.config.email

        pca_evec_file = os.path.join(
            parent_task.params['out_prj_dir'],
            'merged.eigenstrat.pca.evec',
        )

        out_file = os.path.join(basedir, 'merged.eigenstrat.pca.evec.tsv')

        task = {
            'func' : create_evec_data_frame,
            'params' : {
                'in_file' : pca_evec_file,
                'out_file' : out_file,
            },
            'stage_name' : stage,
            'uid' : 'all-chroms',
            'drm_params' :
                to_json(create_evec_data_frame_lsf_params(email)),
            'parents' : [ parent_task ],
        }

        df_task = self.workflow.add_task(**task)

        return df_task

    def create_eigenstrat_smartpca_task(self, parent_task):
        stage = '6-eigenstrat-smartpca'
        basedir = os.path.join(self.config.rootdir, stage)
        email = self.config.email

        ped_file = "{}.ped".format(parent_task.params['out_path'])
        map_file = "{}.map".format(parent_task.params['out_path'])

        task = {
            'func' : eigenstrat_smartpca_analysis,
            'params' : {
                'in_ped_file' : ped_file,
                'in_map_file' : map_file,
                'out_prj_dir' : basedir,
            },
            'stage_name' : stage,
            'uid' : 'all-chroms',
            'drm_params' :
                to_json(eigenstrat_smartpca_analysis_lsf_params(email)),
            'parents' : [ parent_task ],
        }

        eigenstrat_task = self.workflow.add_task(**task)

        return eigenstrat_task

    def create_plink_merge_prune_file_task(self, parent_tasks):
        stage = '5-plink-merge-prune-files'
        basedir = os.path.join(self.config.rootdir, stage)
        email = self.config.email

        parent_tasks_sorted = sorted(parent_tasks, key=lambda t: t.id)

        first_task = parent_tasks_sorted[0]
        remaining_tasks = parent_tasks_sorted[1:]

        merge_list_file = os.path.join(basedir, 'allfiles.txt')
        self._create_merge_list(merge_list_file, remaining_tasks)

        output_path = os.path.join(basedir, 'merged')

        task = {
            'func' : plink_merge_pruned_files,
            'params' : {
                'in_ref' : first_task.params['out_path'],
                'in_merge_file' : merge_list_file,
                'out_path' : output_path,
            },
            'stage_name' : stage,
            'uid' : 'all-chroms',
            'drm_params' :
                to_json(plink_merge_pruned_files_lsf_params(email)),
            'parents' : parent_tasks_sorted,
        }

        merge_task = self.workflow.add_task(**task)

        return merge_task


    def _create_merge_list(self, merge_file, tasks):
        ensure_directory(os.path.dirname(merge_file))
        with open(merge_file, 'w') as f:
            for t in tasks:
                print(t.params['out_path'], file=f)

    def create_plink_extract_prune_tasks(self, parent_tasks):
        tasks = []
        stage = '4-plink-extract-prune'
        basedir = os.path.join(self.config.rootdir, stage)
        email = self.config.email

        for ptask in sorted(parent_tasks, key=lambda t: t.id):
            plink_extract_file = "{}.prune.in".format(ptask.params['out_path'])
            orig_binary_data = ptask.params['in_path']
            chrom = ptask.params['chrom']
            output_path = os.path.join(basedir, chrom, 'c{}.extracted'.format(chrom))

            task = {
                'func' : plink_extract_prune,
                'params' : {
                    'in_path' : orig_binary_data,
                    'in_extract' : plink_extract_file,
                    'out_path' : output_path,
                    'chrom' : chrom,
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' :
                    to_json(plink_extract_prune_lsf_params(email)),
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_plink_ld_prune_tasks(self, parent_tasks):
        tasks = []
        stage = '3-plink-ld-prune'
        basedir = os.path.join(self.config.rootdir, stage)
        email = self.config.email

        for ptask in sorted(parent_tasks, key=lambda t: t.id):
            chrom = ptask.params['chrom']
            output_path = os.path.join(basedir, chrom, 'c{}-pruned'.format(chrom))

            task = {
                'func' : plink_ld_prune,
                'params' : {
                    'in_path' : ptask.params['out_path'],
                    'out_path' : output_path,
                    'chrom' : chrom,
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' :
                    to_json(plink_ld_prune_lsf_params(email)),
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_plink_binary_tasks(self, parent_tasks):
        tasks = []
        stage = '2-plink-binaries'
        basedir = os.path.join(self.config.rootdir, stage)
        email = self.config.email

        for ptask in sorted(parent_tasks, key=lambda t: t.id):
            chrom = ptask.params['chrom']
            output_path = os.path.join(basedir, chrom, 'c{}'.format(chrom))

            task = {
                'func' : plink_binary,
                'params' : {
                    'in_vcf' : ptask.params['out_vcf'],
                    'out_path' : output_path,
                    'chrom' : chrom,
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' :
                    to_json(plink_binary_lsf_params(email)),
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

    def create_filter_biallelic_snps_tasks(self):
        tasks = []
        stage = '1-filter-biallelic-snps'
        basedir = os.path.join(self.config.rootdir, stage)
        email = self.config.email

        for chrom in self.config.chroms:
            vcf = self.config.vcfs[chrom]
            output_vcf = 'filtered.snps.c{chrom}.vcf.gz'.format(chrom=chrom)
            task = {
                'func'   : filter_biallelic_snps,
                'params' : {
                    'chrom' : chrom,
                    'in_vcf' : vcf,
                    'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                    'in_min_vqslod' : self.config.vqslod_threshold,
                },
                'stage_name' : stage,
                'uid' : '{chrom}'.format(chrom=chrom),
                'drm_params' :
                    to_json(filter_biallelic_snps_lsf_params(email)),
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

def create_evec_data_frame(in_file, out_file):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/pca/make-pca-evec-data-frame.py'),
    }

    cmd_args = merge_params(default, args)

    cmd = ( "python -u {script} "
            "--src={in_file} "
            "--out={out_file} ").format(**cmd_args)

    return cmd

def create_evec_data_frame_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "hall-lab",
        'M' : 4000000,
        'R' : 'select[mem>4000] rusage[mem=4000]',
    }

def eigenstrat_smartpca_analysis(in_ped_file, in_map_file, out_prj_dir):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/pca/eigenstrat.mk'),
    }

    cmd_args = merge_params(default, args)

    cmd = ( "make -f {script} "
            "INPUT_PED={in_ped_file} "
            "INPUT_MAP={in_map_file} "
            "PRJ_DIR={out_prj_dir}" ).format(**cmd_args)

    return cmd

def eigenstrat_smartpca_analysis_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "hall-lab",
        'M' : 16000000,
        'R' : 'select[mem>16000] rusage[mem=16000]',
    }

def plink_merge_pruned_files(in_ref, in_merge_file, out_path):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/pca/plink-merge-pruned.sh'),
    }

    cmd_args = merge_params(default, args)

    cmd = ("{script} {in_ref} {in_merge_file} {out_path}").format(**cmd_args)

    return cmd

def plink_merge_pruned_files_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "long",
        'M' : 16000000,
        'R' : 'select[mem>16000] rusage[mem=16000]',
    }

def plink_extract_prune(in_path, in_extract, out_path, **kwargs):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/pca/plink-extract-prune.sh'),
    }

    cmd_args = merge_params(default, args)

    cmd = ("{script} {in_path} {in_extract} {out_path}").format(**cmd_args)

    return cmd

def plink_extract_prune_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "long",
        'M' : 16000000,
        'R' : 'select[mem>16000] rusage[mem=16000]',
    }

def plink_ld_prune(in_path, out_path, **kwargs):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/pca/plink-ld-prune.sh'),
    }

    cmd_args = merge_params(default, args)

    cmd = ("{script} {in_path} {out_path}").format(**cmd_args)

    return cmd

def plink_ld_prune_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "long",
        'M' : 16000000,
        'R' : 'select[mem>16000] rusage[mem=16000]',
    }

def plink_binary(in_vcf, out_path, **kwargs):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/pca/create-plink-binary.sh'),
    }

    cmd_args = merge_params(default, args)

    cmd = ("{script} {in_vcf} {out_path}").format(**cmd_args)

    return cmd

def plink_binary_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "long",
        'M' : 8000000,
        'R' : 'select[mem>8000] rusage[mem=8000]',
    }

def filter_biallelic_snps(in_vcf, out_vcf, in_min_vqslod, **kwargs):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/pca/vcf-filter-biallelic-snps.sh'),
    }

    cmd_args = merge_params(default, args)

    cmd = ("{script} {in_vcf} {out_vcf} {in_min_vqslod}").format(**cmd_args)

    return cmd

def filter_biallelic_snps_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "long",
        'M' : 8000000,
        'R' : 'select[mem>8000] rusage[mem=8000]',
    }
