#!/usr/bin/env python

import sys
import copy
import numpy as np
from scipy.spatial.distance import cdist
from itertools import product

from libpdb import read_pdb, match_pdb, pdb_to_R
from libsup import ls_rmsd

PARAM_SPHERE_RADII=6.0
PARAM_CUTOFF_s=[2.0, 4.0]

def SphereGrinderSingle(pdb_fn, ref, Rref, sphere_s, residue_s,\
        ignore_chain=False, use_calpha=False, cutoff_s=PARAM_CUTOFF_s):
    #
    pdb = read_pdb(pdb_fn, ignore_chain=ignore_chain, use_calpha=use_calpha)
    status = match_pdb(ref, pdb, pdb_fn)
    #
    sg = {}
    sg_global = 0.0
    #
    Rpdb = pdb_to_R(pdb, residue_s=residue_s)
    if status:
        for k,sphere in enumerate(sphere_s):
            rmsd = []
            for s in sphere:
                rmsd.append(ls_rmsd(Rref[s], Rpdb[s]))
            rmsd = min(rmsd)
            sg[residue_s[k]] = rmsd
            for cutoff in cutoff_s:
                if rmsd < cutoff:
                    sg_global += 1
    else:
        ref_tmp = copy.deepcopy(ref)
        Rref_common, sphere_s_common = get_grinder_sphere(ref_tmp, residue_s=residue_s, pdb=pdb)
        for k,sphere in enumerate(sphere_s):
            frac = float(len(sphere[0]))/float(len(sphere_s[k][0]))
            rmsd = []
            for s in sphere:
                rmsd.append(ls_rmsd(Rref_common[s], Rpdb[s]))
            rmsd = min(rmsd) / frac
            sg[residue_s[k]] = rmsd
            for cutoff in cutoff_s:
                if rmsd < cutoff:
                    sg_global += frac
    #
    sg_global /= float(len(residue_s))*len(cutoff_s)/100.0
    return sg_global, sg

def get_grinder_sphere(ref, residue_s, pdb=None):
    Rref = pdb_to_R(ref, residue_s=residue_s, ref=pdb)
    n_atm = [0] + [len(ref[residue]) for residue in residue_s]
    n_atm = np.cumsum(n_atm)
    for i,residue in enumerate(residue_s):
        ref[residue].set_index(n_atm[i])
    #
    sphere_s = []
    for residue in residue_s:
        dist_s = cdist([ref[residue].R("CA")], Rref)
        sphere = np.where(dist_s < PARAM_SPHERE_RADII)[1]
        #
        i_res = np.unique(np.array([np.where(i_atm >= n_atm)[0][-1] for i_atm in sphere]))
        swap = np.array([ref[residue_s[i]].swap for i in i_res])
        i_swap = np.array(np.where(swap)[0])
        if len(i_swap) == 0:
            sphere_s.append([sphere])
        else:
            sphere_tmp = []
            for prod in product([False, True], repeat=len(i_swap)):
                prod = i_res[i_swap[np.array(prod)]]
                if len(prod) == 0:
                    sphere_tmp.append(sphere)
                else:
                    index = []
                    for i,residue in enumerate(residue_s):
                        if i in prod:
                            index.append(ref[residue].swapindex)
                        else:
                            index.append(ref[residue].i_atm)
                    index = np.concatenate(index)
                    sphere_tmp.append(index[sphere])
            sphere_s.append(np.unique(sphere_tmp, axis=0))
    return Rref, sphere_s

def SphereGrinder(ref_fn, pdb_fn_s, ignore_chain=False, use_calpha=False, cutoff_s=PARAM_CUTOFF_s):
    ref = read_pdb(ref_fn, ignore_chain=ignore_chain, use_calpha=use_calpha)
    residue_s = sorted(ref.keys())
    Rref, sphere_s = get_grinder_sphere(ref, residue_s)
    #
    score_s = []
    for pdb_fn in pdb_fn_s:
        try:
            score_s.append(SphereGrinderSingle(pdb_fn, ref, Rref, sphere_s, residue_s,\
                    ignore_chain=ignore_chain, use_calpha=use_calpha, cutoff_s=cutoff_s))
        except:
            score_s.append((0.0, {}))
    return score_s, residue_s
