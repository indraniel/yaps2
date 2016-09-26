import os, pwd
import pkg_resources
from cosmos.api import Cosmos, Dependency, default_get_submit_args
from yaps2.utils import to_json, merge_params, natural_key

class Config(object):
    def __init__(self, job_db, input_vcfs_file,
                       percentiles_file, samples_file, tranches_file, plink_fam_file,
                       project_name, email, workspace):
        self.email = email
        self.db = job_db
        self.project_name = project_name
        self.rootdir = workspace
        self.control_samples_file = samples_file
        self.plink_fam_file = plink_fam_file

        self.ensure_rootdir()

        if self.email is None:
            self.email = self.setup_email()

        if self.db is None:
            self.db = os.path.join(
                os.path.abspath(self.rootdir),
                '.job_queue.db'
            )

        self.collect_tranche_values(tranches_file)
        self.tranche_intervals = self.setup_tranche_intervals()

        self.collect_percentiles(percentiles_file)

        self.vcfs = self.collect_input_vcfs(input_vcfs_file)
        self.chroms = self.get_ordered_chroms()

    def ensure_rootdir(self):
        if not os.path.exists(self.rootdir):
            os.makedirs(self.rootdir)

    def setup_email(self):
        user_id = pwd.getpwuid( os.getuid() ).pw_name
        return '{}@genome.wustl.edu'.format(user_id)

    def collect_tranche_values(self, infile):
        with open(infile, 'r') as f:
            rawdata = [ tuple(line.rstrip().split("\t")) for line in f ]

        self.tranches = {}
        for category in ('snps', 'indels'):
            tranch_levels = { int(x[1]) : float(x[2]) for x in rawdata if x[0] == category }
            self.tranches[category] = tranch_levels

    def setup_tranche_intervals(self):
        tranche_intervals = {}
        for category in ('snps', 'indels') :
            tranches = self.tranches[category]
            tranche_intervals[category] = dict(
                zip(
                    (1,2,3),
                    ( (tranches[1], 100000),
                      (tranches[2], tranches[1]),
                      (tranches[3], tranches[2]) )
                )
            )
        return tranche_intervals

    def collect_percentiles(self, infile):
        # expecting a tsv file of "<snp-or-indel>\t<percentile>\t<vqslod min>\t<vqslod max>" lines
        with open(infile, 'r') as f:
            rawdata = [ tuple(line.rstrip().split("\t")) for line in f ]

        self.percentiles = {}
        for category in ('snps', 'indels') :
            percentiles = { int(x[1]) : ( float(x[2]), float(x[3]) ) for x in rawdata if x[0] == category }
            self.percentiles[category] = percentiles

    def collect_input_vcfs(self, infile):
        # expecting a tsv file of <chrom>\t<path-to-vcf-file> lines
        with open(infile, 'r') as f:
            vcfs = [ tuple(line.rstrip().split("\t")) for line in f ]
        return dict(vcfs)

    def get_ordered_chroms(self):
        chroms = sorted(self.vcfs.keys(), key=natural_key)
        return chroms

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
        partition_tasks = self.create_vcf_partition_tasks()
        plink_pipeline_tasks = self.create_plink_pipeline_tasks(partition_tasks)

    def create_plink_pipeline_tasks(self, parent_tasks):
        tasks = []
        stage = '2-plink-pipeline'
        basedir = os.path.join(self.config.rootdir, stage)
        email = self.config.email

        for ptask in sorted(parent_tasks, key=lambda t: t.id):
            chrom = ptask.params['in_chrom']
            label = ptask.params['in_label']
            method = ptask.params['in_method']
            category = ptask.params['in_type']

            output_dir = os.path.join(basedir, category, method, label, chrom)
#            ensure_directory(output_dir)

            task = {
                'func' : plink_pipeline,
                'params' : {
                    'in_vcf' : ptask.params['out_vcf'],
                    'in_trio_fam' : self.config.plink_fam_file,
                    'chrom' : chrom,
                    'type' : category,
                    'method' : method,
                    'chrom' : chrom,
                    'label' : label,
                    'out_dir' : output_dir,
                },
                'stage_name' : stage,
                'uid' : '{category}:{method}:{label}:{chrom}'.format(
                    chrom=chrom, method=method, category=category, label=label
                ),
                'drm_params' :
                    to_json(plink_pipeline_lsf_params(email)),
                'parents' : [ptask],
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

        pass

    def create_vcf_partition_tasks(self):
        all_tasks = []
        cases = (
            ('tranche', 'snps', self.config.tranche_intervals['snps']),
            ('tranche', 'indels', self.config.tranche_intervals['indels']),
            ('percentile', 'snps', self.config.percentiles['snps']),
            ('percentile', 'indels', self.config.percentiles['indels']),
        )

        for case in cases:
            tasks = self.generate_vcf_partition_tasks(*case)
            all_tasks.extend( tasks )

        return all_tasks

    def generate_vcf_partition_tasks(self, method, category, intervals):
        # method: 'tranche' or 'percentile'
        # category: 'snps' or 'indels'
        # label: 
        #     tranche : 1, 2 or 3
        #     percentile : 10, 20, 30, 40, 50, 60, 70, 80, 90, 100
        tasks = []
        for label in sorted(intervals.keys()) :
            interval = intervals[label]
            partition_tasks = self.create_vcf_partition_chromosome_tasks(
                method=method,
                label=str(label),
                category=category,
                interval=interval,
            )
            tasks.extend(partition_tasks)

        return tasks

    def create_vcf_partition_chromosome_tasks(self, method, label, category, interval):
        tasks = []
        stage = '1-partition-vcfs'
        basedir = os.path.join(self.config.rootdir, stage, category, method, label)
        email = self.config.email

        for chrom in self.config.chroms:
            vcf = self.config.vcfs[chrom]
            output_vcf = 'selected.c{chrom}.vcf.gz'.format(chrom=chrom)
            task = {
                'func'   : vcf_partition,
                'params' : {
                    'in_vcf' : vcf,
                    'out_vcf' : os.path.join(basedir, chrom, output_vcf),
                    'in_min_vqslod' : interval[0],
                    'in_max_vqslod' : interval[1],
                    'in_samples' : self.config.control_samples_file,
                    'in_type': category,
                    'in_method': method,
                    'in_chrom' : chrom,
                    'in_label' : label,
                },
                'stage_name' : stage,
                'uid' : '{category}:{method}:{label}:{chrom}'.format(
                    chrom=chrom, method=method, category=category, label=label
                ),
                'drm_params' :
                    to_json(vcf_partition_lsf_params(email)),
            }
            tasks.append( self.workflow.add_task(**task) )

        return tasks

def vcf_partition(in_vcf, out_vcf, in_min_vqslod, in_max_vqslod, in_samples, in_type, in_method, in_chrom, in_label):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/mie/vcf-partition.sh'),
    }

    cmd_args = merge_params(default, args)

    cmd = ( "{script} {in_vcf} {out_vcf} "
            "{in_min_vqslod} "
            "{in_max_vqslod} "
            "{in_samples} "
            "{in_type} "
            "{in_method}").format(**cmd_args)

    return cmd

def vcf_partition_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "long",
        'M' : 8000000,
        'R' : 'select[mem>8000] rusage[mem=8000]',
    }

def plink_pipeline(in_vcf, in_trio_fam, out_dir, **kwargs):
    args = locals()
    default = {
        'script' : pkg_resources.resource_filename('yaps2', 'resources/mie/plink.mk'),
    }

    cmd_args = merge_params(default, args)

    cmd = ( "make -f {script} "
            "INPUT_VCF={in_vcf} "
            "TRIO_FAM={in_trio_fam} "
            "PRJ_DIR={out_dir}" ).format(**cmd_args)

    return cmd

def plink_pipeline_lsf_params(email):
    return  {
        'u' : email,
        'N' : None,
        'q' : "long",
        'M' : 8000000,
        'R' : 'select[mem>8000] rusage[mem=8000]',
    }
