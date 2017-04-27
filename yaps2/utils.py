import json, re, os

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

class Region(object):
    def __init__(self, reference_index, string):
        self._load_chromosomes(reference_index)
        self.chrom = None
        self.start = None
        self.end = None
        self._parse_region(string)

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

