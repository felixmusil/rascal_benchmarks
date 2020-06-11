from .prettyjson import prettyjson
import json

def tojson(fn, data):
    with open(fn, 'w') as f:
        data_pretty = prettyjson(data,indent=2, maxlinelength=80)
        f.write(data_pretty)

def fromjson(fn):
    with open(fn, 'r') as f:
        data = json.load(f)
    return data