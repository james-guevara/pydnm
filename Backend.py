#!/usr/bin/env python3
import os, sys
def err_fh(f):
    """checks if a file exists"""
    if not os.path.isfile(f):
        print('FATAL ERROR: {} does not exist'.format(f))
        sys.stderr.write('FATAL ERROR: {} does not exist\n'.format(f))
        sys.exit(1)
def get_path():
    try:
        root = __file__
        if os.path.islink(root): root = os.path.realpath(root)
        return os.path.dirname(os.path.abspath(root))
    except:
        sys.stderr.write('FATAL ERROR {} NOT FOUND\n'.format(__file__))
        sys.exit(1)
def pseudoautosome(gen):
    par = {}
    # hg38 0-base positions ^o^
    par['chrX'] = [ [10000,2781479],[155701381,156030895] ]
    par['chrY'] = [ [10000,2781479],[56887902,57217415] ]
    if gen in ['hg19', 'b37','hg37']:
        par['chrX'] = [ [600001,2699520],[154931043,155260560] ]
        par['chrY'] = [ [10000,2649520],[59034049,59363566] ]
    return par
def tokenize(l):
    """splits a line by spaces or tabs, returns a list object"""
    r=l.rstrip().split('\t')
    if len(r)==1: r=l.rstrip().split(' ')
    if len(r)==1:
        sys.stderr.write('WARNING: {} not tokenized correctly. Check if {} is Tab or Space delimited\n'.format(l,l))
        return 0
    else: return r
def try_index(l,ind):
    v=None
    try: v=l[ind]
    except IndexError: v=None
    return v