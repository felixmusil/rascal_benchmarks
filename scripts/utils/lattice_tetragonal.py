#  Adapted from https://github.com/libAtoms/silicon-testing-framework

# This script defines a test case which computes one or more physical
# properties with a given model
#
# INPUTS:
#   model.calculator -- an ase.calculator.Calculator instance
#     this script can assume the calculator is checkpointed.
#
# OUTPUTS:
#   properties -- dictionary of key/value pairs corresponding
#     to physical quantities computed by this test

# standard ASE structure generation routines
from ase.units import GPa

import ase.io, sys, ase.constraints

# set of utility routines specific this this model/testing framework
from .utilities import relax_atoms_cell, evaluate
import matscipy.elasticity

from ase.optimize import FIRE
from ase.optimize.precon.precon import Exp
from ase.optimize.precon.lbfgs import PreconLBFGS
# import model

import numpy as np

class FixedLine_forces_only(ase.constraints.FixedLine):
    def adjust_positions(self, atoms, newpositions):
        pass

    def adjust_forces(self, atoms, forces):
        forces[self.a] = self.dir * np.dot(forces[self.a], self.dir)

def do_lattice(bulk, use_precon=True, elastic=True, tol=1.0e-3):

   print ("unrelaxed bulk")
   ase.io.write(sys.stdout, bulk, format='extxyz')

   # use one of the routines from utilities module to relax the initial
   # unit cell and atomic positions
   if use_precon:
       bulk = relax_atoms_cell(bulk, tol=tol, traj_file="bulk.relax.extxyz", symmetrize=True)
   else:
       bulk = relax_atoms_cell(bulk, tol=tol, traj_file=None, method='fire', symmetrize=True)

   print ("relaxed bulk")
   ase.io.write(sys.stdout, bulk, format='extxyz')

   print ("calculating elastic constants")
   precon = None
   if use_precon:
      precon = Exp(3.0)
   opt = lambda atoms, **kwargs: PreconLBFGS(atoms, precon=precon, **kwargs)
   if elastic:
       # reset calculator to non-symmetrized one (not optimal, but would otherwise need to have optimizer used by fit_elastic_constants to reset symmetry for each relaxation):w
       calc = bulk.get_calculator().calc
       bulk.set_calculator(calc)
       try:
          elastic_consts = matscipy.elasticity.fit_elastic_constants(bulk, symmetry='tetragonal_high', optimizer=opt, logfile=sys.stdout)
       except RuntimeError:
          # fallback on FIRE if we get a linesearch failure with LBFGS
          opt = FIRE
          elastic_consts = matscipy.elasticity.fit_elastic_constants(bulk, symmetry='tetragonal_high', optimizer=opt, logfile=sys.stdout)

       c11 = elastic_consts[0][0,0]/GPa
       c33 = elastic_consts[0][2,2]/GPa
       c12 = elastic_consts[0][0,1]/GPa
       c13 = elastic_consts[0][0,2]/GPa
       c44 = elastic_consts[0][3,3]/GPa
       c66 = elastic_consts[0][5,5]/GPa

   print ("calculating E vs. V")
   V0 = bulk.get_volume()
   dV = bulk.get_volume()*0.025
   E_vs_V=[]
   scaled_bulk = bulk.copy()
   scaled_bulk.set_calculator(bulk.get_calculator())
   print ("bulk going into E vs. V")
   ase.io.write(sys.stdout, scaled_bulk, format='extxyz')
   f = open("relaxed_E_vs_V_configs.xyz", "w")

   scaled_bulk = bulk.copy()
   scaled_bulk.set_calculator(bulk.get_calculator())
   constraints = []

   cell_scalings = np.linspace(0.9**(1.0/3.0), 1.1**(1.0/3.0), 30)
   for i,cell_scaling in enumerate(cell_scalings):
      scaled_bulk = bulk.copy()
      scaled_bulk.set_calculator(bulk.get_calculator())
      print ("doing volume step",i)
      scaled_bulk.set_cell(scaled_bulk.get_cell()*cell_scaling, scale_atoms=True)

      scaled_bulk = relax_atoms_cell(scaled_bulk, tol=tol, traj_file=None, constant_volume=True, method='fire', symmetrize=True, max_steps=500)

      print ("done relaxing step",i)
      ase.io.write(f, scaled_bulk, format='extxyz')
      f.flush()
      E_vs_V.insert(0, (scaled_bulk.get_volume()/len(bulk), scaled_bulk.get_potential_energy()/len(bulk)) )
      evaluate(scaled_bulk)
      print ("done evaluate step",i)
      print ("EV ",i, scaled_bulk.get_volume(), scaled_bulk.get_potential_energy(), scaled_bulk.get_stress())

   for (V, E) in E_vs_V:
     print ("EV_final ", V, E)

   if elastic:
       return (c11, c33, c12, c13, c44, c66, E_vs_V)
   else:
       return (E_vs_V)
