#!/usr/bin/env python

import numpy as np
from scipy.spatial.distance import cdist

from libpdb import read_pdb, match_pdb, pdb_to_R
from libsup import ls_rmsd

PARAM_SPHERE_RADII=6.0
PARAM_CUTOFF_s=[2.0, 4.0]

def SphereGrinderSingle(pdb_fn, ref, Rref, residue_s,\
        ignore_chain=False, use_calpha=False):
    pdb = read_pdb(pdb_fn, ignore_chain=ignore_chain, use_calpha=use_calpha)
    status = match_pdb(ref, pdb, pdb_fn)
    #
    sg = {}
    sg_global = 0.0
    #
    Rpdb = pdb_to_R(pdb, residue_s=residue_s)
    if status:
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
    else:
        Rref_common = pdb_to_R(ref, residue_s=residue_s, ref=pdb)
        for residue in residue_s:
            dist_s = cdist([ref[residue].R("CA")], Rref_common)
            sphere = np.where(dist_s < PARAM_SPHERE_RADII)[1]
            Sref = Rref_common[sphere]
            Spdb = Rpdb[sphere]
            #
            rmsd = ls_rmsd(Spdb, Sref)
            sg[residue] = rmsd
            for cutoff in PARAM_CUTOFF_s:
                if rmsd < cutoff:
                    sg_global += 1
    #
    sg_global /= float(len(residue_s))*len(PARAM_CUTOFF_s)/100.0
    return sg_global, sg

def SphereGrinder(ref_fn, pdb_fn_s, ignore_chain=False, use_calpha=False):
    ref = read_pdb(ref_fn, ignore_chain=ignore_chain, use_calpha=use_calpha)
    residue_s = sorted(ref.keys())
    Rref = pdb_to_R(ref, residue_s=residue_s)
    #
    score_s = []
    for pdb_fn in pdb_fn_s:
        score_s.append(SphereGrinderSingle(pdb_fn, ref, Rref, residue_s,\
                ignore_chain=ignore_chain, use_calpha=use_calpha))
    return score_s, residue_s
