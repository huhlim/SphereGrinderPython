#!/usr/bin/env python

import sys
import numpy as np

GREEK_s = ['A','B','G','D','E','Z','H']
BACKBONE_s = ['N','CA','C','O']
OXT_s = ['OXT','O1','O2','OT1','OT2']

def ATOM_SORT_RULE(atm_1, atm_2):
    if atm_1[0] in BACKBONE_s:
        ia_1 = BACKBONE_s.index(atm_1[0])
    else:
        ia_1 = 4
    if atm_2[0] in BACKBONE_s:
        ia_2 = BACKBONE_s.index(atm_2[0])
    else:
        ia_2 = 4
    if ia_1 != ia_2:    # backbone atoms
        return cmp(ia_1, ia_2)
    #
    ia_1 = GREEK_s.index(atm_1[0][1])
    ia_2 = GREEK_s.index(atm_2[0][1])
    #
    if len(atm_1[0]) > 2:
        ia_1 += 0.1*float(atm_1[0][2])
    if len(atm_2[0]) > 2:
        ia_2 += 0.1*float(atm_2[0][2])
    #
    return cmp(ia_1, ia_2)

class Residue:
    def __init__(self, resNo, chain):
        self.resNo = resNo
        self.chain = chain
        self.atom_s = []
    def __repr__(self):
        return '%s %s %2d'%(self.resNo, self.chain, len(self.atom_s))
    def put_ATOM(self, atmName, R):
        self.atom_s.append((atmName, R))
    def sort(self):
        self.atom_s.sort(ATOM_SORT_RULE)
    def resort(self, atmName_s):
        atom_s = []
        atmName_curr = self.atmName()
        for atmName in atmName_s:
            atom_s.append(self.atom_s[atmName_curr.index(atmName)])
        self.atom_s = atom_s
        self._atmName = atmName_s
    def atmName(self):
        if '_atmName' not in dir(self):
            self._atmName = []
            for atom in self.atom_s:
                self._atmName.append(atom[0])
        return self._atmName
    def R(self, atmName=None):
        if atmName != None:
            for atom in self.atom_s:
                if atom[0] == atmName:
                    return atom[1]
        elif '_R' in dir(self):
            return self._R
        else:
            self._R = []
            for atom in self.atom_s:
                self._R.append(atom[1])
            self._R = np.array(self._R)
            return self._R

def read_pdb(pdb_fn):
    pdb = {}
    with open(pdb_fn) as fp:
        for line in fp:
            if not line.startswith("ATOM"):
                continue
            atmName = line[13:16].strip()
            if atmName in OXT_s : continue
            if atmName[0] == 'H': continue
            #
            resNo = line[22:27]
            chain = line[21].replace(" ",'_')
            #
            key = (resNo, chain)
            if key not in pdb:
                pdb[key] = Residue(resNo, chain)
            #
            R = np.array([line[30:38], line[38:46], line[46:54]], dtype=float)
            pdb[key].put_ATOM(atmName, R)
    for residue in pdb:
        pdb[residue].sort()
    return pdb

def match_pdb(ref, pdb, pdb_fn):
    for residue in ref:
        if residue not in pdb:
            sys.stderr.write("ERROR: residue %s %s is not exists in the PDB file %s\n"%(residue[0], residue[1], pdb_fn))
            return False
        if not set(ref[residue].atmName()).issubset(set(pdb[residue].atmName())):
            sys.stderr.write("ERROR: residue %s %s is not same in the PDB file %s\n"%(residue[0], residue[1], pdb_fn))
            return False
        if ref[residue].atmName() == pdb[residue].atmName():
            continue
        pdb[residue].resort(ref[residue].atmName())
    return True

def pdb_to_R(pdb, residue_s=[]):
    if len(residue_s) == 0:
        residue_s = list(pdb.keys())
    R = []
    for residue in residue_s:
        R.append(pdb[residue].R())
    return np.concatenate(R)
