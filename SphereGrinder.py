#!/usr/bin/env python

import sys
import argparse
from libsg import SphereGrinder

def main():
    arg = argparse.ArgumentParser(description='SphereGrinder')
    arg.add_argument('-r', '--ref', dest='ref_fn', metavar='REFERENCE',\
            help='reference protein structure', required=True)
    arg.add_argument('-m', '--model', dest='pdb_fn_s', metavar='MODEL',\
            help='model protein structure', required=True, nargs='+')
    arg.add_argument('-f', '--full', dest='full', \
            help='print full log', default=False, action='store_true')
    arg.add_argument('--ignore-chain', dest='ignore_chain',\
            help='ignore chain IDs', default=False, action='store_true')
    #
    if len(sys.argv) == 1:
        return arg.print_help()
    arg = arg.parse_args()
    #
    score_s, residue_s =\
            SphereGrinder(arg.ref_fn, arg.pdb_fn_s, ignore_chain=arg.ignore_chain)
    #
    sys.stdout.write("# reference: %s\n"%arg.ref_fn)
    sys.stdout.write("#\n")
    #
    for i,pdb_fn in enumerate(arg.pdb_fn_s):
        sys.stdout.write("SphereGrinderGlobal %6.2f %s\n"%(score_s[i][0], pdb_fn))
        if arg.full:
            for residue in residue_s:
                sys.stdout.write("SphereGrinderLocal  %7.3f %s %s %s\n"%\
                        (score_s[i][1][residue], residue[0], residue[1], pdb_fn))
            sys.stdout.write("#\n")

if __name__=='__main__':
    main()
