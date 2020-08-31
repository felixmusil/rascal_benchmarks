
import numpy as np
import sys
from os.path import join,dirname


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

from concurrent.futures import as_completed, ProcessPoolExecutor, ThreadPoolExecutor


def init_proc(sparse_points_fn, zeta):
    global X_pseudo
    global soap_nograd,soap_grad
    global kernel_grad,kernel_nograd
    # compute_gradient = True
    # zeta = 4
    # sparse_points_fn = sp_fn
    X_pseudo = load_obj(sparse_points_fn)
    hypers = X_pseudo.representation._get_init_params()
    hypers['compute_gradients'] = False
    soap_nograd = SphericalInvariants(**hypers)
    hypers['compute_gradients'] = True
    soap_grad = SphericalInvariants(**hypers)
    kernel_grad = Kernel(soap_grad, name='GAP', zeta=zeta, target_type='Structure', kernel_type='Sparse')
    kernel_nograd = Kernel(soap_nograd, name='GAP', zeta=zeta, target_type='Structure', kernel_type='Sparse')

def compute_gap(frame, compute_gradient, compute_stress, zeta):
    if compute_gradient:
        feat = soap_grad.transform([frame])
        en_row = kernel_grad(feat, X_pseudo)
    else:
        feat = soap_nograd.transform([frame])
        en_row = kernel_nograd(feat, X_pseudo)
    grad_rows = None
    virial_rows = None
    # if compute_gradient and compute_stress:
    #     grad_rows = kernel_grad(feat, X_pseudo, grad=(True, False), compute_stress=True)
    #     virial_rows = grad_rows[-6:]
    #     grad_rows = grad_rows[:-6]
    if compute_gradient and not compute_stress:
        # grad_rows = kernel_grad(feat, X_pseudo, grad=(True, False), compute_stress=False)
        grad_rows = kernel_grad(feat, X_pseudo, grad=(True, False))
    feat = []
    return en_row, grad_rows, virial_rows


class KnmPool(object):
    def __init__(self, ncpu=1, energy_tag=None, forces_tag=None, stress_tag=None):
        self.ncpu = ncpu
        self.energy_tag = energy_tag
        self.forces_tag = forces_tag
        self.stress_tag = stress_tag

        if self.energy_tag is None:
            raise ValueError('energy_tag should be provided')

    def prepare_run(self,frames,self_contributions,n_sparse,n_frames,n_stress):
        self.Y0 = np.zeros((n_frames,1))
        self.n_atoms = np.zeros((n_frames,1))
        self.n_atoms_v = []
        self.compute_gradient = []
        self.compute_virial = []
        self.energies = []
        self.forces = []
        self.stress = []
        self.n_grad_stride = [0]
        self.n_virial_stride = [0]
        self.frame_ids = {'energy':[],'forces':[],'stress':[]}
        for i_frame, frame in enumerate(frames):
            numbers = frame.get_atomic_numbers()
            self.n_atoms[i_frame] = len(frame)

            for sp in numbers:
                self.Y0[i_frame] += self_contributions[sp]

            if self.energy_tag not in frame.info:
                raise ValueError(
                    'Could not find: "{}" in frame {}'.format(self.energy_tag, i_frame))

            self.energies.append(frame.info[self.energy_tag])
            self.frame_ids['energy'].append(i_frame)
            if self.forces_tag in frame.arrays:
                self.compute_gradient.append(True)
                self.forces.extend(frame.get_array(self.forces_tag).flatten()[:, None])
                self.n_grad_stride.append(len(frame)*3)
                self.frame_ids['forces'].append(i_frame)
                if self.stress_tag in frame.info:
                    self.compute_virial.append(True)
                    stress = -frame.get_volume()* full_3x3_to_voigt_6_stress(frame.info[self.stress_tag])
                    self.stress.extend(stress.reshape((n_stress, 1)))
                    self.n_virial_stride.append(n_stress)
                    self.n_atoms_v.extend([len(frame)]*n_stress)
                    self.frame_ids['stress'].append(i_frame)
                else:
                    self.compute_virial.append(False)
                    self.n_virial_stride.append(0)

            else:
                self.compute_gradient.append(False)
                self.n_grad_stride.append(0)

        self.n_grad_stride = np.cumsum(self.n_grad_stride)
        self.n_virial_stride = np.cumsum(self.n_virial_stride)

        self.energies = np.asarray(self.energies)[:,None]-self.Y0
        self.grads = -np.asarray(self.forces)
        self.m_virials = -np.asarray(self.stress)
        self.n_atoms_v = np.asarray(self.n_atoms_v).reshape((-1,1))

    def run(self, frames, zeta, sparse_points_fn, self_contributions):
        X_pseudo = load_obj(sparse_points_fn)
        n_sparse = X_pseudo.size()
        n_frames = len(frames)
        n_stress = 6
        self.prepare_run(frames,self_contributions,n_sparse,n_frames,n_stress)

        KNM_e = np.ones((n_frames, n_sparse))
        if self.forces_tag is not None:
            KNM_f = np.ones((self.n_grad_stride[-1], n_sparse))
        if self.stress_tag is not None:
            KNM_v = np.ones((self.n_virial_stride[-1], n_sparse))
        # compute_gradient = True
        with ThreadPoolExecutor(max_workers=self.ncpu, initializer=init_proc,
                                 initargs=(sparse_points_fn, zeta)) as executor:
            future_to_compute = {executor.submit(
                        compute_gap, frame, compute_gradient, compute_virial, zeta) : i_frame
                        for i_frame,(frame, compute_gradient, compute_virial) in enumerate(zip(
                                                    frames, self.compute_gradient, self.compute_virial))}
            pbar = tqdm(total=len(future_to_compute),leave=False, desc='KNM')
            for future in as_completed(future_to_compute):
                i_frame = future_to_compute[future]
                en_row, grad_rows, virial_rows = future.result()
                KNM_e[i_frame] = en_row
                if grad_rows is not None:
                    KNM_f[self.n_grad_stride[i_frame]:self.n_grad_stride[i_frame+1]] = grad_rows
                if virial_rows is not None:
                    KNM_v[self.n_virial_stride[i_frame]:self.n_virial_stride[i_frame+1]] = virial_rows

                pbar.update()
        result = {'energy':{'KNM':KNM_e, 'y':self.energies, 'n_atoms':self.n_atoms}}
        if self.forces_tag is not None:
            result['grads'] = {'KNM':KNM_f, 'y':self.grads}
        if self.stress_tag is not None:
            result['stress'] = {'KNM':KNM_v, 'y':self.m_virials, 'n_atoms':self.n_atoms_v}
        return result

    def srun(self, frames, zeta, sparse_points_fn, self_contributions):
        init_proc(sparse_points_fn, zeta)
        global X_pseudo
        X_pseudo = load_obj(sparse_points_fn)
        n_sparse = X_pseudo.size()
        n_frames = len(frames)
        n_stress = 6
        self.prepare_run(frames,self_contributions,n_sparse,n_frames,n_stress)

        KNM_e = np.ones((n_frames, n_sparse))
        if self.forces_tag is not None:
            KNM_f = np.ones((self.n_grad_stride[-1], n_sparse))
        if self.stress_tag is not None:
            KNM_v = np.ones((self.n_virial_stride[-1], n_sparse))
        pbar = tqdm(total=n_frames,leave=False, desc='KNM')
        for i_frame, frame in enumerate(frames):
            en_row, grad_rows, virial_rows = compute_gap(frame, self.compute_gradient[i_frame],
                                                         self.compute_virial[i_frame], zeta)
            KNM_e[i_frame] = en_row
            if grad_rows is not None:
                KNM_f[self.n_grad_stride[i_frame]:self.n_grad_stride[i_frame+1]] = grad_rows
            if virial_rows is not None:
                KNM_v[self.n_virial_stride[i_frame]:self.n_virial_stride[i_frame+1]] = virial_rows
            pbar.update()
        result = {'energy':{'KNM':KNM_e, 'y':self.energies, 'n_atoms':self.n_atoms}}
        if self.forces_tag is not None:
            result['grads'] = {'KNM':KNM_f, 'y':self.grads}
        if self.stress_tag is not None:
            result['stress'] = {'KNM':KNM_v, 'y':self.m_virials, 'n_atoms':self.n_atoms_v}
        return result

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
            kernel = self.kernel(managers, self.X_train, (compute_gradients, False), compute_stress)
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

def train_gap_model(kernel, X_pseudo, energy, self_contributions, grads=None, stress=None,
                    lambdas=None, jitter=1e-8):
    KMM = kernel(X_pseudo)
    Y_e = energy['y'].reshape((-1, 1))
    KNM_e = energy['KNM']
    Natoms = energy['n_atoms']

    delta = np.std(Y_e)
    # lambdas[0] is provided per atom hence the '* np.sqrt(Natoms)'
    # the first n_centers rows of KNM are expected to refer to the
    # property
    reg = 1 / (lambdas[0] / delta * np.sqrt(Natoms))
    Y = np.dot(KNM_e.T, Y_e * reg**2)
    K = KMM + np.dot(KNM_e.T , KNM_e* reg**2)

    if grads is not None:
        reg = 1 / (lambdas[1] / delta)
        Y_f = grads['y'].reshape((-1, 1))
        KNM_f = grads['KNM']
        Y += np.dot(KNM_f.T, Y_f* reg**2)
        K += np.dot(KNM_f.T, KNM_f* reg**2)
    if stress is not None:
        n_atoms = stress['n_atoms'].copy()
        reg = 1 / (lambdas[2] / delta * np.sqrt(n_atoms))
        Y_v = stress['y'].reshape((-1, 1))
        KNM_v = stress['KNM']
        Y += np.dot(KNM_v.T, Y_v* reg**2)
        K += np.dot(KNM_v.T, KNM_v* reg**2)

    eig,_ = np.linalg.eig(K)
    print('jitter ',jitter, eig.min())
    K[np.diag_indices_from(K)] += jitter

    # weights = np.linalg.solve(K, Y)
    weights,_,_,_ = np.linalg.lstsq(K, Y, rcond=None)
    model = KRR(weights, kernel, X_pseudo, self_contributions)

    return model