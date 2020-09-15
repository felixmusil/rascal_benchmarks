
import numpy as np
import sys
from os.path import join,dirname
import ase

sys.path.insert(0, join(dirname(__file__), '../'))
from path import STRUCTURE_PATH, RASCAL_BUILD_PATH, BUILD_PATH

sys.path.insert(0, RASCAL_BUILD_PATH)
from rascal.representations import SphericalInvariants, SphericalExpansion
from rascal.representations.spherical_invariants import get_power_spectrum_index_mapping
from rascal.representations import SphericalInvariants
from rascal.models import Kernel, SparsePoints
from rascal.neighbourlist import AtomsList
from rascal.utils import from_dict, CURFilter, FPSFilter, fps, to_dict
from rascal.utils.random_filter import RandomFilter
from rascal.utils.io import dump_obj,load_obj,BaseIO

from tqdm import tqdm


from .io import fromfile, tofile


class KRR(BaseIO):
    """Kernel Ridge Regression model. Only compatible fully with sparse GPR
    training for the moment.

    Parameters
    ----------
    weights : np.array
        weights of the model

    kernel : Kernel
        kernel class used to train the model

    X_train : PseudoPoints
        reference samples used for the training

    self_contributions : dictionary
        map atomic number to the property baseline, e.g. isolated atoms
        energies when the model has been trained on total energies.
    """

    def __init__(self, weights, kernel, X_train, self_contributions):
        super(KRR, self).__init__()
        # Weights of the krr model
        self.weights = weights
        self.kernel = kernel
        self.X_train = X_train
        self.self_contributions = self_contributions
        self.target_type = kernel.target_type

    def _get_property_baseline(self, managers):
        """build total baseline contribution for each prediction"""
        if self.target_type == 'Structure':
            Y0 = np.zeros(len(managers))
            for i_manager, manager in enumerate(managers):
                if isinstance(manager, ase.Atoms):
                    numbers = manager.get_atomic_numbers()
                    for sp in numbers:
                        Y0[i_manager] += self.self_contributions[sp]
                else:
                    for at in manager:
                        Y0[i_manager] += self.self_contributions[at.atom_type]
        elif self.target_type == 'Atom':
            n_centers = 0
            for manager in managers:
                n_centers += len(manager)
            Y0 = np.zeros(n_centers)
            i_center = 0
            for manager in managers:
                for sp in manager.get_atomic_numbers():
                    Y0[i_center] = self.self_contributions[sp]
                    i_center += 1
        return Y0

    def _preprocess_input(self, managers, KNM, compute_gradients=False, compute_stress=False):
        """compute prediction kernel and total baseline contributions"""
        from rascal.utils.io import is_npy

        if KNM is not None: # if the KNM matrix is provided
            kernel = KNM
        else: # if the representation is provided
            kernel = self.kernel(managers, self.X_train, (compute_gradients, False),)
        Y0 = self._get_property_baseline(managers)
        return kernel, Y0

    def predict(self, managers, KNM=None, compute_gradients=False, compute_stress=False):
        """Predict properties associated with the atomic structures in managers
        or their derivative w.r.t. atomic positions (if compute_gradients==True).

        Parameters
        ----------
        managers : AtomsList
            list of atomic structures with already computed features compatible
            with representation in kernel
        compute_gradients : bool, optional
            predict the gradients of the property w.r.t atomic positions,
            by default False

        Returns
        -------
        np.array
            predictions
        """
        KNM, Y0 = self._preprocess_input(managers, KNM, compute_gradients, compute_stress)
        if compute_gradients is False:
            return Y0 + np.dot(KNM, self.weights).reshape(-1)
        else:
            aa = np.dot(KNM, self.weights)
            return np.dot(KNM, self.weights)

    def get_weigths(self):
        return self.weights

    def _get_init_params(self):
        init_params = dict(weights=self.weights, kernel=self.kernel,
                           X_train=self.X_train, self_contributions=self.self_contributions)
        return init_params

    def _set_data(self, data):
        pass

    def _get_data(self):
        return dict()

def extract_ref(frames,info_key='energy',array_key='forces'):
    y,f = [], []
    for frame in frames:
        y.append(frame.info[info_key])
        if array_key is None:
            pass
        elif array_key == 'zeros':
            f.append(np.zeros(frame.get_positions().shape))
        else:
            f.append(frame.get_array(array_key))
    y= np.array(y)
    try:
        f = np.concatenate(f)
    except:
        pass
    return y,f

def get_grad_strides(frames):
    Nstructures = len(frames)
    Ngrad_stride = [0]
    Ngrads = 0
    n_atoms = []
    for frame in frames:
        n_at = len(frame)
        n_atoms.append(n_at)
        Ngrad_stride.append(n_at*3)
        Ngrads += n_at*3
    Ngrad_stride = np.cumsum(Ngrad_stride)
    return np.array(n_atoms).reshape((-1,1)),Ngrad_stride

def get_strides(frames):
    Nstructures = len(frames)
    Ngrad_stride = [0]
    Ngrads = 0
    for frame in frames:
        n_at = len(frame)
        Ngrad_stride.append(n_at*3)
        Ngrads += n_at*3
    Ngrad_stride = np.cumsum(Ngrad_stride) + Nstructures
    return Nstructures,Ngrads,Ngrad_stride

def _split_KNM(ids, frames, KNM, grad_strides):
    Nstruct, Ngrad, _ = get_strides(frames)
    KNM_ = np.zeros((Nstruct+Ngrad, KNM.shape[1]))
    i_struct,i_grad = 0, Nstruct
    for idx in ids:
        KNM_[i_struct] = KNM[idx]
        kk = KNM[grad_strides[idx]:grad_strides[idx+1]]
        KNM_[i_grad:i_grad+kk.shape[0]] = kk
        i_struct += 1
        i_grad += kk.shape[0]
    return KNM_

def split_KNM(frames,KNM,train_ids, test_ids):
    Nstruct, Ngrad, grad_strides = get_strides(frames)

    frames_train = [frames[ii] for ii in train_ids]
    frames_test = [frames[ii] for ii in test_ids]
    y_train, f_train = extract_ref(frames_train,'dft_energy','dft_force')
    y_test, f_test = extract_ref(frames_test,'dft_energy','dft_force')

    KNM_train = _split_KNM(train_ids, frames_train, KNM, grad_strides)
    KNM_test = _split_KNM(test_ids, frames_test, KNM, grad_strides)
    return (frames_train, y_train, f_train, KNM_train), (frames_test, y_test, f_test, KNM_test)


def train_gap_model_0(kernel, frames, KNM_, X_pseudo, y_train, self_contributions, grad_train=None,
                    lambdas=None, jitter=1e-8):
    KMM = kernel(X_pseudo)
    Y = y_train.reshape((-1, 1)).copy()
    KNM = KNM_.copy()
    n_centers = Y.shape[0]
    Natoms = np.zeros(n_centers)
    Y0 = np.zeros((n_centers, 1))
    for iframe, frame in enumerate(frames):
        Natoms[iframe] = len(frame)
        numbers = frame.get_atomic_numbers()
        for sp in numbers:
            Y0[iframe] += self_contributions[sp]
    Y = Y - Y0
    delta = np.std(Y)

    KNM[:n_centers] /= lambdas[0] / delta * np.sqrt(Natoms)[:, None]
    Y /= lambdas[0] / delta * np.sqrt(Natoms)[:, None]

    if grad_train is not None:
        KNM[n_centers:] /= lambdas[1] / delta
        F = grad_train.reshape((-1, 1)).copy()
        F /= lambdas[1] / delta
        Y = np.vstack([Y, F])

    K = KMM + np.dot(KNM.T, KNM)
    eig,_ = np.linalg.eig(K)
    if eig.min() < 0:
        jitter = 1.05*abs(eig.min())

    # print('jitter ',jitter, eig.min())
    K[np.diag_indices_from(K)] += jitter

    Y = np.dot(KNM.T, Y)
    weights,_,_,_ = np.linalg.lstsq(K, Y, rcond=None)
    model = KRR(weights, kernel, X_pseudo, self_contributions)

    # avoid memory clogging
    del K, KMM
    K, KMM = [], []

    return model


def get_cv_scores(cv,KNM, frames, kernel, X_pseudo,self_contributions,lamda_es,lamda_fs,jitter=1e-7, **kwargs):
    Natoms = []
    y_baseline = []
    for frame in frames:
        Natoms.append(len(frame))
        y_baseline.append(len(frame)*self_contributions[14])
    Natoms = np.array(Natoms).reshape(-1)
    y_baseline = np.array(y_baseline).reshape(-1)
    n_atoms, grad_strides = get_grad_strides(frames)
    scores = []
    for train, val in tqdm(cv.split(np.ones((len(frames),1))), leave=False, total=cv.n_splits):
        (frames_train, y_train, f_train, KNM_train), (frames_test,
                                    y_test, f_test, KNM_test) = split_KNM(frames,KNM, train, val)
        for lamda_e in tqdm(lamda_es,leave=False):
            for lamda_f in tqdm(lamda_fs,leave=False):
                if lamda_e>lamda_f: continue
                model = train_gap_model_0(kernel, frames_train, KNM_train, X_pseudo, y_train, self_contributions,
                        grad_train=-f_train, lambdas=[lamda_e,lamda_f], jitter=jitter)
                frames_v = [frames[ids] for ids in val]
                yp = model.predict(frames_v, KNM=KNM_test[:len(frames_test)])
                fp = -model.predict(frames_v, KNM=KNM_test[len(frames_test):], compute_gradients=True)

                en_score = get_score((yp.flatten()-y_baseline[val])/Natoms[val],
                                     (y_test.flatten()-y_baseline[val])/Natoms[val])
                score = {k+'_e':v for k,v in en_score.items()}
                f_score = get_score(fp.flatten(), f_test.flatten())
                score.update({k+'_f':v for k,v in f_score.items()})
                score.update(lambda_e=lamda_e,lambda_f=lamda_f, **kwargs)
                scores.append(score)
    return scores