#!/usr/bin/env python

import numpy as np
import copy

class Quaternion:
    def __init__(self, q):
        self.q = q
    def __repr__(self):
        return self.q
    def rotate(self):
        if 'R' in dir(self):
            return self.R
        #
        self.R = np.zeros((3,3))
        #
        self.R[0][0] = self.q[0]**2 + self.q[1]**2 - self.q[2]**2 - self.q[3]**2
        self.R[0][1] = 2.0*(self.q[1]*self.q[2] - self.q[0]*self.q[3])
        self.R[0][2] = 2.0*(self.q[1]*self.q[3] + self.q[0]*self.q[2])
        #
        self.R[1][0] = 2.0*(self.q[1]*self.q[2] + self.q[0]*self.q[3])
        self.R[1][1] = self.q[0]**2 - self.q[1]**2 + self.q[2]**2 - self.q[3]**2
        self.R[1][2] = 2.0*(self.q[2]*self.q[3] - self.q[0]*self.q[1])
        #
        self.R[2][0] = 2.0*(self.q[1]*self.q[3] - self.q[0]*self.q[2])
        self.R[2][1] = 2.0*(self.q[2]*self.q[3] + self.q[0]*self.q[1])
        self.R[2][2] = self.q[0]**2 - self.q[1]**2 - self.q[2]**2 + self.q[3]**2
        return self.R

def v_size(v):
    return sqrt(v.dot(v))

def v_norm(v):
    return v/v_size(v)

def ls_rmsd(_X, _Y):    # _X -> model ; _Y -> ref
    X = copy.copy(_X)
    Y = copy.copy(_Y)
    n = float(len(X))
    #
    X_cntr = X.transpose().sum(1)/n
    Y_cntr = Y.transpose().sum(1)/n
    X -= X_cntr
    Y -= Y_cntr
    Xtr = X.transpose()
    Ytr = Y.transpose()
    X_norm = (Xtr*Xtr).sum()
    Y_norm = (Ytr*Ytr).sum()
    #
    Rmatrix = np.zeros(9).reshape((3,3))
    for i in range(3):
        for j in range(3):
            Rmatrix[i][j] = Xtr[i].dot(Ytr[j])
    #
    S = np.zeros(16).reshape((4,4))
    S[0][0] =  Rmatrix[0][0] + Rmatrix[1][1] + Rmatrix[2][2]
    S[1][0] =  Rmatrix[1][2] - Rmatrix[2][1]
    S[0][1] =  S[1][0]
    S[1][1] =  Rmatrix[0][0] - Rmatrix[1][1] - Rmatrix[2][2]
    S[2][0] =  Rmatrix[2][0] - Rmatrix[0][2]
    S[0][2] =  S[2][0]
    S[2][1] =  Rmatrix[0][1] + Rmatrix[1][0]
    S[1][2] =  S[2][1]
    S[2][2] = -Rmatrix[0][0] + Rmatrix[1][1] - Rmatrix[2][2]
    S[3][0] =  Rmatrix[0][1] - Rmatrix[1][0]
    S[0][3] =  S[3][0]
    S[3][1] =  Rmatrix[0][2] + Rmatrix[2][0]
    S[1][3] =  S[3][1]
    S[3][2] =  Rmatrix[1][2] + Rmatrix[2][1]
    S[2][3] =  S[3][2]
    S[3][3] = -Rmatrix[0][0] - Rmatrix[1][1] + Rmatrix[2][2] 
    #
    eigl,eigv = np.linalg.eigh(S)
    q = eigv.transpose()[-1]
    sU = Quaternion(q).rotate()
    sT = Y_cntr - sU.dot(X_cntr)
    #
    rmsd = np.sqrt(max(0.0, (X_norm + Y_norm - 2.0*eigl[-1]))/n)
    return rmsd
    ##
    #X = copy.copy(_X)
    #for i in range(len(X)):
    #    X[i] = sU.dot(X[i]) + sT
    #dXY = np.sqrt(((_Y-X)**2).sum(1))
    #return rmsd, dXY
