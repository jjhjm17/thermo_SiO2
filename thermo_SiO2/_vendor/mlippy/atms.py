# This file is derived from MLIP-2 (mlippy).
#
# Copyright (c) 2023 Alexander Shapeev (Skoltech)
#
# no modification
#
# See LICENSE.mlippy for the original BSD 2-Clause license.
from __future__ import print_function
import ase
# from cfgs import *
from .cfgs import *

def ase_atom(cfg):
    atom = ase.Atoms(positions=cfg.pos,cell=cfg.lat,numbers=cfg.types,
                    pbc=[True,True,True])
    if cfg.energy is not None:
        atom.energy = cfg.energy
    if cfg.forces is not None:
        atom.forces = cfg.forces
    if cfg.stresses is not None:
        atom.stresses = cfg.stresses
    return atom

def ase_cfg(atom):
    cfg = Cfg()
    cfg.pos = atom.positions
    cfg.lat = atom.cell
    cfg.types = atom.numbers
    if hasattr(atom, 'energy'):
        cfg.energy = atom.energy
    if hasattr(atom, 'forces'):
        cfg.forces = atom.forces
    if hasattr(atom, 'stresses'):    
        cfg.stresses = atom.stresses
    return cfg
    
class ase_atomparser:
    def __init__(self, file, max_cfgs = None):
        self.atoms = []
        self.file = file
        self.max_cfgs = max_cfgs
    def __enter__(self):
        while True:
            if self.max_cfgs is not None and len(self.atoms) == self.max_cfgs:
                break
            cfg = readcfg(self.file)
            if cfg.types is not None :
                atom = ase_atom(cfg)
                self.atoms.append(atom)
            else:
                break
        return self.atoms
    def __exit__(self, *args):
        self.cfgs = []


def ase_loadcfgs(filename,max_cfgs = None):
    with open(filename, 'r') as file:
        with ase_atomparser(file,max_cfgs) as atoms:
            return atoms
        
def ase_savecfgs(filename,atoms,desc = None):
    with open(filename, 'w') as file:
        for atom in atoms:
            cfg = ase_cfg(atom)
            savecfg(file,cfg,desc)
            print ("",file=file)

