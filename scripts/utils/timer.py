from timeit import default_timer as timer
from contextlib import contextmanager
import numpy as np

class Timer(object):
    def __init__(self, tag='', logger=None):
        self.tag = tag
        self.elapsed = []
        self.start = None
        self.end = None
    def __enter__(self):
        self.start = timer()
    def __exit__(self, type, value, traceback):
        self.end = timer()
        self.elapsed.append(self.end-self.start)
        
    def mean(self):
        return np.mean(self.elapsed)
    def stdev(self):
        return np.std(self.elapsed)
    def min(self):
        return np.min(self.elapsed)
    def max(self):
        return np.max(self.elapsed)
    def samples(self):
        return self.elapsed
    def dumps(self):
        data = dict(tag=self.tag,mean=self.mean(),stdev=self.stdev(),
                    min=self.min(), max=self.max(),
                    samples=self.samples())
        return data