# pyDNM
----

Classifier used to discover single-nucleotide variant (SNV) de novo mutations (DNMs) and indel DNMs.

---

:warning: Still under construction :construction:
USE AT YOUR OWN RISK

---

## Install

```
$ pip install https://github.com/dantaki/pydnm/releases/download/v0.1.0.0/pyDNM-0.1.0.0.tar.gz 
```

---

## Usage

```
usage: 
                        ooooooooo   oooo   oooo  oooo     oooo 
oooooooooo oooo   oooo   888    88o  8888o  88    8888o   888  
 888    888 888   888    888    888  88 888o88    88 888o8 88  
 888    888  888 888     888    888  88   8888    88  888  88  
 888ooo88      8888     o888ooo88   o88o    88   o88o  8  o88o 
 888           888
o888        o8o888                                         

Python port of forestDNM for SNVs+INDELs : http://sebatlab.ucsd.edu/software-data 
Version 0.1.0.0    Authors: Danny Antaki, Aojie Lian, James Guevara    
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

```


