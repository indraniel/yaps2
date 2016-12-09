from __future__ import print_function, division

import os, sys, pwd, json
import pkg_resources
from cosmos.api import Cosmos, Dependency, default_get_submit_args
from yaps2.utils import to_json, merge_params, natural_key, ensure_directory

class Config(object):
    def __init__(self, job_db,
                 input_json_sample_bams, project_name, email, workspace, drm_job_group):
        self.email = email
        self.db = job_db
        self.project_name = project_name
        self.rootdir = workspace
        self.input_json_sample_bams = input_json_sample_bams
        self.drm_job_group = drm_job_group

        self.ensure_rootdir()

        if self.email is None:
            self.email = self.setup_email()

        if self.db is None:
            self.db = os.path.join(
                os.path.abspath(self.rootdir),
                '.job_queue.db'
            )

        self.sample_data = self.collect_sample_data()

    def ensure_rootdir(self):
        if not os.path.exists(self.rootdir):
            os.makedirs(self.rootdir)

    def setup_email(self):
        user_id = pwd.getpwuid( os.getuid() ).pw_name
        return '{}@genome.wustl.edu'.format(user_id)

    def collect_sample_data(self):
        # expecting the JSON to have this structure
        # {
        #    "<organism-sample-id>": {
        #        "administration-project": "CCDG Testing - Gold Whole Genomes with TruSeq", 
        #        "bams": [
        #            "/path/to/gerald/bam1.bam", 
        #            "/path/to/geral/bam2.bam"
        #        ], 
        #        "meta": {
        #            "ethnicity": null, 
        #            "gender": "male", 
        #            "nomenclature": "CCDG", 
        #            "original-name": "H_XYZ-sample-1234"
        #        }
        #    }, 
        # }
        #
        # NOTE: the 'adminstration-project' key can be substituted by 'analysis-project'
        with open(self.input_json_sample_bams, 'r') as json_data:
            d = json.load(json_data)
        return d

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
        speedseq_tasks = self.create_speedseq_realign_tasks()

    def create_speedseq_realign_tasks(self):
        tasks = []
        stage = '1-exec-speedseq-realign'
        basedir = os.path.join(self.config.rootdir, stage)
        email = self.config.email
        lsf_job_group = self.config.drm_job_group
        sample_data = self.config.sample_data

        for sample_id in sample_data.keys():
            bam_paths = sample_data[sample_id]['bams']
            sample_name = sample_data[sample_id]['meta']['original-name']
            output_prefix = os.path.join(basedir, sample_id, "{}.b38.realign".format(sample_id))
            tmpdir = os.path.join(basedir, sample_id, 'tmpdir')
            input_bams = ' '.join(bam_paths)

            task = {
                'func'   : exec_speedseq,
                'params' : {
                    'output_prefix' : output_prefix,
                    'tmpdir' : tmpdir,
                    'input_bams' : input_bams,
                },
                'stage_name' : stage,
                'uid' : sample_id,
                'drm_params' :
                    to_json(exec_speedseq_lsf_params(email, lsf_job_group)),
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

def exec_speedseq(output_prefix, tmpdir, input_bams, **kwargs):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/b38/speedseq-realign.sh'),
        # build 38
        'reference' : os.path.join(
            '/gscmnt/gc2802/halllab/ccdg_resources/genomes',
            'human/GRCh38DH/bwa/0_7_12/all_sequences.fa'
        ),
    }

    cmd_args = merge_params(default, args)

    cmd = ("{script} {output_prefix} {tmpdir} {reference} {input_bams}").format(**cmd_args)

    return cmd

def exec_speedseq_lsf_params(email, job_group):
    return  {
        'g' : job_group,
        'u' : email,
        'N' : None,
        'q' : "long",
        'M' : 50000000, # 50_000_000 (50 GB)
        'R' : 'select[mem>45000] rusage[mem=48000] span[hosts=1]',
        'n' : 8
    }
