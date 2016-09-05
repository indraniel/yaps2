import json, re

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
