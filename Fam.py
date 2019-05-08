#!/usr/bin/env python3
from pyDNM.Backend import err_fh,tokenize
import sys
class Fam():
    def __init__(self):
        self.offspring = {} # [child ID] = (Father,Mother)
        self.sex = {} # [child ID] = Sex as a string!
    def load_fam(self,fh=None):
        err_fh(fh)
        with open(fh) as f:
            for l in f:
                r = tokenize(l)
                if r==0: continue
                if len(r)<5:
                    sys.stderr.write('WARNING: {} does not contain 5 elements\n'.format(l))
                    continue
                iid,sex = str(r[1]),str(r[4])
                dad,mom = str(r[2]),str(r[3])
                if sex!='1' and sex!='2':
                    sys.stderr.write('WARNING: {} not accepted format for sex. Error found in {}\n'.format(sex,iid))
                    continue
                self.sex[iid]=sex
                if dad == '0' or mom =='0': continue
                # process offspring only
                self.offspring[iid]=(dad,mom)
