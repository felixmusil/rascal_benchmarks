import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def make_colormap(seq,name='CustomMap',ncolors=256):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap(name, cdict,N=ncolors)


def get_energyCmap(Ncolors=256):
    mycolors = {'orange':np.array((255,128,0))/255.,'purple':np.array((88,24,112))/255.,'deepblue':np.array((75,53,182))/255.,
                                'blue':np.array((88,126,227))/255.,'lightblue':np.array((155,209,252))/255.,'white':np.array((245,245,245))/255.}

    seq = [ mycolors['orange'],mycolors['purple'],0.25,mycolors['purple'],mycolors['deepblue'],0.5,mycolors['deepblue'],mycolors['blue'],
           0.75,mycolors['blue'],mycolors['lightblue']]

    return make_colormap(seq,name='MCenergy',ncolors=Ncolors)

def get_periodicCmap(Ncolors=256):
    mycolors = {'white':np.array((213,200,209))/255.,'red':np.array((228,77,12))/255.,'middle':np.array((64,44,57))/255.,
                                                                                                  'blue':np.array((12,98,206))/255.}
    seq = [ mycolors['white'],mycolors['red'],0.25,mycolors['red'],mycolors['middle'],0.5,mycolors['middle'],mycolors['blue'],
       0.75,mycolors['blue'],mycolors['white']]

    return make_colormap(seq,name='MCperio',ncolors=Ncolors)