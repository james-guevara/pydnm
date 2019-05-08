#!/usr/bin/env python3
from pyDNM.Backend import err_fh,tokenize,try_index
from pyDNM.Features import Feature
import numpy as np
import sys

def genotype_pl_index(gt=None):
    """
        In presence of the GT field the same
        ploidy is expected and the canonical order is used; without GT field, diploidy is assumed.  If A is the allele in
        REF and B,C,...  are the alleles as ordered in ALT, the ordering of genotypes for the likelihoods is given by:
        
        F(j/k) = (k*(k+1)/2)+j 
    """
    _a = gt.replace('|','/').split('/')
    if len(_a)==2:
        a1,a2 = int(_a[0]),int(_a[1])
        t = a1
        if a2 < a1: a1,a2=a2,t
        return int((a2*(a2+1)/2.)+a1)
    else:
        return -1

def intersect_range(x,y):
    # x and y are lists with two elements
    # the 0-based start and end; x = [99, 100] for example 
    return len(list(range(max(x[0],y[0]), min(x[1],y[1]))))

class Vcf():
    def __init__(self):
        self.ids={} # [iid] = column index
        self.format=None # [format] = index
        self.gt=None # [iid] = genotype
        self.missing=None# True if one sample has a missing genotype
    def index_samples(self,l=None,Fam=None):
        """
        Processes CHROM header line
        stores column indices with sample ID
        """
        r = tokenize(l)
        if r==0:
            print('FATAL ERROR: unable to tokenize header {}\n'.format(l))
            sys.stderr.write('FATAL ERROR: unable to tokenize header {}\n'.format(l))
            sys.exit(1)
        for i in range(9,len(r)):
            if Fam.sex.get(r[i])==None: continue
            self.ids[r[i]]=i
    def index_format(self,l=None):
        """
        Determines indices of FORMAT values
        """
        self.format={}
        r = l.split(':')
        for i in range(0,len(r)): self.format[str(r[i])]=i
    def load_genotypes(self,r=None):
        """
        Stores genotypes
        """
        self.missing = [False, []]
        self.gt={}
        for iid in self.ids:
            self.gt[iid]=str(r[self.ids[iid]].split(':').pop(self.format['GT']))
        for x in self.gt:
            if '.' in self.gt[x]:
                self.missing[0]=True
                self.missing[1].append('{}:{}'.format(x,self.gt[x]))
    def check_genotypes(self,iid=None,variant=None):
        """
        Logic to check for genotype value in dict
        """
        gt=None
        if self.gt.get(iid)!=None: gt=self.gt[iid]
        if gt==None:
            sys.stderr.write('WARNING: missing genotype entry: {} {}\n'.format(variant,iid))
        return gt
    def allele_depth(self,entry=None,dnm=None,par=None):
        """
        returns 
            the allele ratio of the dnm and parent alleles
            the log2 coverage ratio (relative to median)
        """
        _buff = 1 # add one to avoid 0 values: taken from forestDNM
        _ad = str(entry).split(':')
        if len(_ad)<self.format['AD']: return -1,np.nan
        else:
            _ad = _ad[self.format['AD']].split(',')
            _dnm = try_index(_ad,dnm)
            _par = try_index(_ad,par)
            if _dnm==None or _par==None: return -1,np.nan
            else: 
                _par = float(_par)+_buff
                ar = float(_dnm)/ (_par)
                med = np.median([float(x) for x in _ad if x!='.']) 
                if not np.isfinite(med): return -1, np.nan
                else: return ar, np.log2( (float(_dnm)+float(_par)) / (med+_buff) )
    def genotype_quals(self,entry=None):
        """
        returns GQ values
        """
        _gq = str(entry).split(':')
        if len(_gq)<self.format['GQ']: return -1
        else:
            _gq = _gq[self.format['GQ']]
            if _gq == '.': return -1
            else: return float(_gq) 
    def phred_quals(self,entry=None,gt=None):
        """
        returns PL values
        """
        _pl = str(entry).split(':')
        if len(_pl)<self.format['PL']: return -1
        else:
            _pl = _pl[self.format['PL']].split(',')
            idx = genotype_pl_index(gt)
            if idx==-1: return -1
            else:
                pl = try_index(_pl,idx)
                if pl==None: return -1
                elif pl=='.': return -1
                else: return float(pl)

    #---------------------------------------------------------------
    # Main Function to parse VCF and extract genotypes
    def parse(self,fh=None, Fam=None,verb=None,ofh=None,pseud=None):
        err_fh(fh)
        vcf_fh=None
        if fh.endswith('.gz'):
            import gzip
            vcf_fh = gzip.open(fh,'rt', 9)
        else: vcf_fh = open(fh,'r')
        if vcf_fh==None:
            print('FATAL ERROR: {} file format unknown.\nAccepted VCF formats: text, bgzip (ends with \'.gz\')\n'.format(fh))
            sys.stderr.write('FATAL ERROR: {} file format unknown.\nAccepted VCF formats: text, bgzip (ends with \'.gz\')\n'.format(fh))
            sys.exit(1)
        out = open(ofh,'w')
        out.write('chrom\tpos\tid\tref\talt\tiid\toffspring_gt\tfather_gt\tmother_gt\t')
        out.write('{}\n'.format(Feature().header()))
        for l in vcf_fh:
            if l.startswith('#CHROM'): self.index_samples(l,Fam)
            if l.startswith('#'): continue
            r=tokenize(l)
            if r==0: continue
            variant = tuple(map(str,r[0:5]))
            # load genotypes
            self.index_format(r[8])
            self.load_genotypes(r)
            chrom_tmp = variant[0]
            if variant[0] != "chrX" and variant[0] != "X" and variant[0] != "chrY" and variant[0] != "Y":
                if self.missing[0]==True:
                    if verb==True: sys.stderr.write('WARNING: missing genotypes {} {}\n'.format(variant,','.join(self.missing[1])))
                    continue
            """
            Foreach trio
            """
            for kid in Fam.offspring:

                # skip if child is female and variant is on Y chromosome
                if Fam.sex[kid] == '2' and r[0].endswith('Y'): continue

                if self.ids.get(kid)==None: continue
                dad,mom = Fam.offspring[kid]
                kgt,dgt,mgt = self.check_genotypes(kid,variant),self.check_genotypes(dad,variant),self.check_genotypes(mom,variant)
                # skip if genotypes are not available
                if kgt==None or dgt==None or mgt==None: continue
                alleles = kgt.replace('|','/').split('/')
                
                # determine which of the offspring alleles are de novo and inherited
                
                # 0-base positions
                pos0 = [int(r[1])-1,int(r[1])+len(r[3])-1]
                chrom = r[0]
                if not r[0].startswith('chr'): chrom = 'chr'+chrom
                par1,par2= None,None #pseudoautosome regions
                if chrom == 'chrX' or chrom=='chrY':
                    par1,par2 = pseud[chrom]
                is_diploid=True
                # skip if male, on sex chromosome, not in PAR, and het genotype
                if (Fam.sex[kid] == '1' and 
                    (chrom == 'chrX' or chrom=='chrY') and 
                    intersect_range(pos0,par1)==0 and 
                    intersect_range(pos0,par2)==0):

                        is_diploid=False # haploid genotype
                if is_diploid==False and len(set(alleles)) > 1: continue

                par,dnm = None,None # parent allele, de novo allele
                for a in alleles:
                    # Diploid
                        # if the offspring allele is NOT in parents genotype ==> De Novo
                    if is_diploid==True:
                        if "." in a: continue
                        if str(a) not in dgt and str(a) not in mgt: dnm=int(a)
                        # if the offspring allele is in ONE of the parents genotype ==> Inherited
                        if str(a) in dgt or str(a) in mgt: par=int(a)
                    # Haploid
                    else:
                        if "." in a: continue
                        if chrom == 'chrX' and str(a) not in mgt: dnm=int(a)
                        if chrom == 'chrY' and str(a) not in dgt: dnm=int(a)
                # For male sex chromosome variants, get the correct parent allele
                if is_diploid==False and dnm != None:
                    if chrom == "chrX": par_alleles = mgt.replace('|','/').split('/')
                    elif chrom == 'chrY': par_alleles = dgt.replace('|','/').split('/')
                    if "." in par_alleles: continue
                    for a in par_alleles:
                        if str(a) not in kgt: par=int(a)

                # skip if there is not one inherited and one de novo allele
                if par==None or dnm==None: continue
                # init Features
                Feat = Feature()
                # load INFO features
                Feat.parse(r)
                # skip multiallelic
                Feat.n_alt = len(str(r[4]).split(','))
                if Feat.n_alt > 1: continue
                #-------------------------------------------------------------------------------------------
                # Allele depth
                if self.format.get('AD')==None: 
                    if verb==True: sys.stderr.write('WARNING: missing allele depth AD {}\n'.format(variant))
                    continue
                #print(l)
                #print(kid)
                if is_diploid==False:
                    if chrom.endswith("X"):
                        dad_ad = -1
                        dad_dp = -1
                        mom_ad,mom_dp = self.allele_depth(r[self.ids[mom]],dnm,par)
                        kid_ad,kid_dp = self.allele_depth(r[self.ids[kid]],dnm,par)
                    else:
                        mom_ad = -1
                        mom_dp = -1
                        dad_ad,dad_dp = self.allele_depth(r[self.ids[dad]],dnm,par)
                        kid_ad,kid_dp = self.allele_depth(r[self.ids[kid]],dnm,par)
                else:
                    if "." in kgt or "." in dgt or "." in mgt: continue
                    dad_ad,dad_dp = self.allele_depth(r[self.ids[dad]],dnm,par)
                    mom_ad,mom_dp = self.allele_depth(r[self.ids[mom]],dnm,par)
                    kid_ad,kid_dp = self.allele_depth(r[self.ids[kid]],dnm,par)
                if is_diploid==True:
                #if variant[0] != "chrX" and variant[0] != "X" and variant[0] != "chrY" and variant[0] != "Y":
                    if dad_ad==-1 or mom_ad==-1 or kid_ad==-1: continue
                # if variant[0] != "chrX" and variant[0] != "X" and variant[0] != "chrY" and variant[0] != "Y":
                    if not np.isfinite(dad_dp) or not np.isfinite(mom_dp) or not np.isfinite(kid_dp): continue

                Feat.p_ar_max= max(dad_ad,mom_ad)
                Feat.p_ar_min= min(dad_ad,mom_ad)
                Feat.o_ar = kid_ad
                Feat.p_dp_max = max(dad_dp,mom_dp)
                Feat.p_dp_min = min(dad_dp,mom_dp)
                Feat.o_dp = kid_dp
                #-------------------------------------------------------------------------------------------
                #-------------------------------------------------------------------------------------------
                # Genotype quality scores
                if self.format.get('GQ')==None:
                    if verb==True: sys.stderr.write('WARNING: missing genotype quality GQ {}\n'.format(variant))
                    continue
                if is_diploid==False:
                    if chrom.endswith("X"):
                        kid_gq,mom_gq = self.genotype_quals(r[self.ids[kid]]),self.genotype_quals(r[self.ids[mom]])
                        dad_gq = -1
                    else:
                        kid_gq,dad_gq = self.genotype_quals(r[self.ids[kid]]),self.genotype_quals(r[self.ids[dad]])
                        mom_gq = -1 
                else:
                    kid_gq,dad_gq,mom_gq = self.genotype_quals(r[self.ids[kid]]),self.genotype_quals(r[self.ids[dad]]),self.genotype_quals(r[self.ids[mom]])
                Feat.p_gq_max = max(dad_gq,mom_gq)
                Feat.p_gq_min = min(dad_gq,mom_gq)
                Feat.o_gq = kid_gq
                #-------------------------------------------------------------------------------------------
                #-------------------------------------------------------------------------------------------
                # Phred scaled genotype likelihoods
                if self.format.get('PL')==None:
                    if verb==True: sys.stderr.write('WARNING: missing Phred-adjusted genotype likelihoods PL {}\n'.format(variant))
                    continue
                if is_diploid==False:
                    if chrom.endswith("X"):
                        kid_pl,mom_pl = self.phred_quals(r[self.ids[kid]],kgt),self.phred_quals(r[self.ids[mom]],mgt)
                        dad_pl = -1
                        kid_m_pl,mom_o_pl = self.phred_quals(r[self.ids[kid]],mgt),self.phred_quals(r[self.ids[mom]],kgt)
                        kid_d_pl = -1
                        dad_o_pl = -1
                    else:
                        kid_pl,dad_pl = self.phred_quals(r[self.ids[kid]],kgt),self.phred_quals(r[self.ids[dad]],dgt)
                        mom_pl = -1
                        kid_d_pl,dad_o_pl = self.phred_quals(r[self.ids[kid]],dgt),self.phred_quals(r[self.ids[dad]],kgt)
                        kid_m_pl = -1
                        mom_o_pl = -1
                else:
                    kid_pl, dad_pl, mom_pl = self.phred_quals(r[self.ids[kid]],kgt),self.phred_quals(r[self.ids[dad]],dgt),self.phred_quals(r[self.ids[mom]],mgt)
                    if kid_pl==-1 or dad_pl==-1 or mom_pl==-1: continue
                    kid_d_pl,kid_m_pl, dad_o_pl, mom_o_pl = self.phred_quals(r[self.ids[kid]],dgt),self.phred_quals(r[self.ids[kid]],mgt),self.phred_quals(r[self.ids[dad]],kgt),self.phred_quals(r[self.ids[mom]],kgt)
                if is_diploid==True:
                    if kid_d_pl==-1 or kid_m_pl==-1 or dad_o_pl==-1 or mom_o_pl == -1: continue
                Feat.p_og_max = max(dad_o_pl,mom_o_pl)
                Feat.p_og_min = min(dad_o_pl,mom_o_pl)
                Feat.p_pg_max = max(dad_pl,mom_pl)
                Feat.p_pg_min = min(dad_pl,mom_pl)
                Feat.og = kid_pl
                if is_diploid==False:
                    if chrom.endswith("X"): Feat.o_pg = kid_m_pl
                    else: Feat.o_pg = kid_d_pl
                else:
                    Feat.o_pg = np.median([kid_d_pl,kid_m_pl])
                #-------------------------------------------------------------------------------------------

                #-------------------------------------------------------------------------------------------
                # Sex Chromosome Methods for Males
                if is_diploid==False:
                        if chrom.endswith('X'):
                            Feat.p_ar_max, Feat.p_ar_min = mom_ad,mom_ad
                            Feat.p_dp_max, Feat.p_dp_min = mom_dp,mom_dp
                            Feat.p_gq_max, Feat.p_gq_min = mom_gq,mom_gq
                            Feat.p_og_max, Feat.p_og_min = mom_o_pl, mom_o_pl
                            Feat.p_pg_max, Feat.p_pg_min = mom_pl, mom_pl
                        else:
                            Feat.p_ar_max, Feat.p_ar_min = dad_ad,dad_ad
                            Feat.p_dp_max, Feat.p_dp_min = dad_dp,dad_dp
                            Feat.p_gq_max, Feat.p_gq_min = dad_gq,dad_gq
                            Feat.p_og_max, Feat.p_og_min = dad_o_pl, dad_o_pl
                            Feat.p_pg_max, Feat.p_pg_min = dad_pl, dad_pl
                #-------------------------------------------------------------------------------------------
                # Output the data
                o = Feat.output()
                out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('\t'.join(variant),kid,kgt,dgt,mgt,o))
        out.close()
