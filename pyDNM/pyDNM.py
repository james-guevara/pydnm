#!/usr/bin/env python3
'''
Copyright <2018> <Danny Antaki>
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''
__version__='0.1.0.0'
from pyDNM.Fam import Fam
from pyDNM.Vcf import Vcf
from pyDNM.Backend import pseudoautosome
from pyDNM.Clf import classify
from argparse import RawTextHelpFormatter
import argparse,os,shutil,sys
__usage__="""
                        ooooooooo   oooo   oooo  oooo     oooo 
oooooooooo oooo   oooo   888    88o  8888o  88    8888o   888  
 888    888 888   888    888    888  88 888o88    88 888o8 88  
 888    888  888 888     888    888  88   8888    88  888  88  
 888ooo88      8888     o888ooo88   o88o    88   o88o  8  o88o 
 888           888
o888        o8o888                                         

python port of forestDNM for SNVs+INDELs : http://sebatlab.ucsd.edu/software-data 
Version {}    Authors: Danny Antaki, Aojie Lian, James Guevara    
                   Contact: dantaki at ucsd dot edu
---------------------------------------------------------------------------------
    pyDNM  -f <in.fam>  -v  <in.vcf>  [-oLgkVh]
    
input arguments:
  
  -v, -vcf    PATH    vcf file
  -f, -fam    PATH    plink fam file
  
optional arguments:

  -o, -out    PATH    output file
  -L, -log    PATH    log file for standard error messages [default: STDOUT]
  -g  -gen    STR     human reference genome version [default: hg38]
  -k                  keep false positive de novos in output
  -V                  print extra warnings
  
  -h, -help           show this message and exit
     
""".format(__version__)
def main():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, usage=__usage__, add_help=False)
    in_args, opt_args = parser.add_argument_group('input arguments'), parser.add_argument_group('optional arguments')
    in_args.add_argument('-v', '-vcf', type=str, default=None, required=True)
    in_args.add_argument('-f', '-fam', type=str, default=None, required=True)
    opt_args.add_argument('-L', '-log', default=None, required=False)
    opt_args.add_argument('-o', '-out', required=False, default="pydnm.out", type=str)
    opt_args.add_argument('-g', '-gen', required=False, default="hg38", type=str,choices=['hg19', 'b37', 'hg37','hg38'])
    opt_args.add_argument('-k', required=False,action="store_true", default=False)
    opt_args.add_argument('-V', required=False, action="store_true", default=False)
    opt_args.add_argument('-h', '-help', required=False, action="store_true", default=False)
    args = parser.parse_args()
    vcf,fam,gen,keep_fp = args.v,args.f,args.g,args.k
    logfh,ofh,verb = args.L,args.o,args.V
    _help = args.h
    if (_help==True or len(sys.argv)==1):
        print(__usage__)
        sys.exit(0)
    if fam==None:
        print('\nFATAL ERROR: fam file required\n')
        sys.exit(1)
    if vcf==None:
        print('\nFATAL ERROR: vcf file required\n')
        sys.exit(1)
    if logfh!=None:
        lfh=open(logfh,'w')
        sys.stderr = lfh
    _Fam = Fam()
    _Fam.load_fam(fam)
    Vcf().parse(vcf,_Fam,verb,ofh,pseudoautosome(gen))
    classify(ofh,keep_fp)
