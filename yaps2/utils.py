import json, re, os, gzip

def to_json(var):
    return json.dumps(var)

def merge_params(common, unique):
    z = common.copy()
    z.update(unique)
    return z

# From: http://stackoverflow.com/questions/2545532/python-analog-of-natsort-function-sort-a-list-using-a-natural-order-algorithm
def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

def ensure_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def empty_gzipped_vcf(path):
    # From http://stackoverflow.com/questions/37874936/how-to-check-empty-gzip-file-in-python
    with gzip.GzipFile(path, 'rb') as f:
        for line in f:
            if not line.startswith('#'):
                return False
    return True

def get_chrom_number(region):
    fmt_chrom = ''
    if ':' in region:
        fmt_chrom = region.split(':')[0]
    else:
        fmt_chrom = region

    chrom = ''
    if 'chr' in fmt_chrom:
        chrom = fmt_chrom.lstrip('chr')
    else:
        chrom = fmt_chrom

    return chrom

class Region(object):
    def __init__(self, reference_index, string):
        self._load_chromosomes(reference_index)
        self.chrom = None
        self.start = None
        self.end = None
        self._parse_region(string)

    @staticmethod
    def _start_lt(start, other_start):
        return start != other_start and (start is None or (other_start is not None and int(start) < int(other_start)))

    @staticmethod
    def _end_lt(end, other_end):
        return end != other_end and (other_end is None or (end is not None and int(end) < int(other_end)))

    def __hash__(self):
        return hash((self.chrom, self.start, self.end))

    def __eq__(self, other):
        return self.chrom == other.chrom and self.start == other.start and self.end == other.end

    def __lt__(self, other):
        return natural_key(self.chrom) < natural_key(other.chrom) or (self.chrom == other.chrom and (self._start_lt(self.start, other.start) or (self.start == other.start and self._end_lt(self.end, other.end))))

    def _load_chromosomes(self, reference_index):
        with open(reference_index, 'r') as fai:
            self.chromosomes = set([ x.strip().split('\t')[0] for x in list(fai)])

    def _parse_region(self, string):
        if string in self.chromosomes:
            # Don't want to accidentally split on colons in chromosome names
            self.chrom = string
        else:
            long_match = re.match('(\S+):(\d+)(?:-(\d+))*', string)
            if long_match:
                self.chrom, self.start, self.end = long_match.group(1, 2, 3)
            else:
                raise RuntimeError('Invalid range: {}'.format(string))

