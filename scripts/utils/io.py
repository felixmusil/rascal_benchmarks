from .prettyjson import prettyjson
import json
import ubjson
import os
import ase
import numpy as np
import pickle

def _decode(o):
    # JSON does not have integer keys so they are converted to string
    # to load the object as it was in python this hook converts to 'int' all
    # dictionary keys that can be converted
    if isinstance(o, str):
        try:
            return int(o)
        except ValueError:
            return o
    elif isinstance(o, dict):
        return {_decode(k): _decode(v) for k, v in o.items()}
    else:
        return o

def topickle(fn, data):
    with open(fn, 'wb') as f:
        pickle.dump(data, f, protocol=4)

def frompickle(fn):
    with open(fn, 'rb') as f:
        data = pickle.load(f)
    return data
    
def tojson(fn, data):
    with open(fn, 'w') as f:
        data_pretty = prettyjson(data,indent=2, maxlinelength=80)
        f.write(data_pretty)

def fromjson(fn):
    with open(fn, 'r') as f:
        data = json.load(f)
    return data

def tofile(fn,frames):
    import numpy as np
    keys = ['positions','cell','numbers','pbc']
    data = {}
    for ii,frame in zip(range(1, len(frames)+1), frames):
        aa = dict(positions=frame.get_positions().tolist(),cell=frame.get_cell().tolist(),
                numbers=frame.get_atomic_numbers().tolist(),pbc=frame.get_pbc().tolist())
        aa['info'] = {}
        for k,v in  frame.info.items():
            if isinstance(v, np.integer):
                aa['info'][k] = int(v)
            elif isinstance(v, np.bool_):
                aa['info'][k] = bool(v)
            elif isinstance(v, np.float):
                aa['info'][k] = float(v)
            elif hasattr(v, 'tolist'):
                aa['info'][k] = v.tolist()
            else:
                aa['info'][k] = v
        aa['arrays'] = {}
        for k,v in frame.arrays.items():
            if k not in keys:
                aa['arrays'][k] = v.tolist()
        data[str(ii)] = aa
    data['ids'] = list(range(1, len(frames)+1))
    data['nextid'] = len(frames)+1
    _, extension = os.path.splitext(fn)
    if extension == '.json':
        with open(fn, 'w') as f:
            # json.dump(data, f)
            data_pretty = prettyjson(data,indent=2, maxlinelength=80)
            f.write(data_pretty)
    elif extension == '.ubjson':
        with open(fn, 'wb') as f:
            ubjson.dump(data, f,no_float32=False)

def fromfile(fn):
    _, extension = os.path.splitext(fn)
    if extension == '.json':
        with open(fn, 'r') as f:
            data = json.load(f)
    elif extension == '.ubjson':
        with open(fn, 'rb') as f:
            data = ubjson.load(f)

    frames = []
    for idx in data['ids']:
        ff = data['{}'.format(idx)]
        frame = ase.Atoms(positions=ff['positions'], cell=ff['cell'],
                            numbers=ff['numbers'], pbc=ff['pbc'])
        if 'info' in ff:
            frame.info = ff['info']
        if 'arrays' in ff:
            for k,v in ff['arrays'].items():
                frame.set_array(k, np.array(v))
        frames.append(frame)
    return frames