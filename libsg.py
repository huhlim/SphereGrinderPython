#!/usr/bin/env python

import numpy as np
from scipy.spatial.distance import cdist

from libpdb import read_pdb, match_pdb, pdb_to_R
from libsup import ls_rmsd

PARAM_SPHERE_RADII=6.0
PARAM_CUTOFF_s=[2.0, 4.0]

def test():
    ref = read_pdb("../test/TR866-D1.pdb")
    pdb = read_pdb("../test/TR866-D1_init.pdb")
    match_pdb(ref, pdb)
    #
    Rref = pdb_to_R(ref)
    Rpdb = pdb_to_R(pdb, residue_s=list(ref.keys()))
    #
    for residue in ref:
        dist_s = cdist([ref[residue].R("CA")], Rref)
        sphere = np.where(dist_s < PARAM_SPHERE_RADII)[1]
        Sref = Rref[sphere]
        Spdb = Rpdb[sphere]
        rmsd,dist = ls_rmsd(Spdb, Sref)
        n_atom = len(sphere)
        n_cut  = [len(np.where(dist < PARAM_CUTOFF)[0]) for PARAM_CUTOFF in PARAM_CUTOFF_s]
        sg = [n_cut[i]/float(n_atom) for i in range(len(PARAM_CUTOFF_s))]
        print(residue, '%6.4f  '%(np.mean(sg)), '%6.4f %6.4f'%tuple(sg))

def SphereGrinderSingle(pdb_fn, ref, Rref, residue_s):
    pdb = read_pdb(pdb_fn)
    status = match_pdb(ref, pdb, pdb_fn)
    if not status:
        sg = {}
        for residue in residue_s:
            sg[residue] = 0.0
        sg_global = 0.0
        return sg_global, sg
    #
    Rpdb = pdb_to_R(pdb, residue_s=residue_s)
    #
    sg = {}
    sg_global = 0.0
    for residue in residue_s:
        dist_s = cdist([ref[residue].R("CA")], Rref)
        sphere = np.where(dist_s < PARAM_SPHERE_RADII)[1]
        Sref = Rref[sphere]
        Spdb = Rpdb[sphere]
        #
        rmsd = ls_rmsd(Spdb, Sref)
        sg[residue] = rmsd
        for cutoff in PARAM_CUTOFF_s:
            if rmsd < cutoff:
                sg_global += 1
    sg_global /= float(len(residue_s))*len(PARAM_CUTOFF_s)/100.0
    return sg_global, sg

def SphereGrinder(ref_fn, pdb_fn_s):
    ref = read_pdb(ref_fn)
    residue_s = sorted(ref.keys())
    Rref = pdb_to_R(ref, residue_s=residue_s)
    #
    score_s = []
    for pdb_fn in pdb_fn_s:
        score_s.append(SphereGrinderSingle(pdb_fn, ref, Rref, residue_s))
    return score_s, residue_s
