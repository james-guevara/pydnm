#!/usr/bin/env python3
import numpy as np
class Feature():
    def __init__(self):
        """
        Below are features taken from the source code of forestDNM
        file R/pDNM.R
        """
        self.filt=None # FILTER
        self.qual=None # QUAL
        # allele ratios
        self.p_ar_max=None # max parent allele ratio (dnm/par)
        self.p_ar_min=None # min parent allele ratio (dnm/par)
        self.o_ar=None # offspring allele ratio (dnm/par)
        # depth
        self.p_dp_max=None # max parent depth
        self.p_dp_min=None # min parent depth
        self.o_dp=None # offspring depth
        # phred lik
        self.p_og_max=None # parent max genotype phred lik for offspring GT
        self.p_og_min=None # parent min genotype phred lik for offspring GT
        self.p_pg_max = None  # parent max genotype phred lik for parent GT
        self.p_pg_min = None  # parent min genotype phred lik for parent GT
        self.og=None # offspring genotype phred lik
        self.o_pg=None # offspring genotype phred lik for parent GT
        # gq
        self.p_gq_max = None  # max parent GQ
        self.p_gq_min = None  # min parent GQ
        self.o_gq = None  # offspring GQ
        self.n_alt = None  # number of alt alleles
        # INFO features
        self.info_feats=[]
        self.info_keys=['VQSLOD','ClippingRankSum','BaseQRankSum','FS','SOR','MQ','MQRankSum','QD','ReadPosRankSum']

    def info_features(self,l=None):
        info={}
        for x in l.split(';'):
            if '=' not in x: continue            
            key,val = x.split('=')
            info[str(key)]=val
        for key in self.info_keys:
            if info.get(key)!=None: self.info_feats.append(float(info[key]))
            else: self.info_feats.append(np.nan)

    def parse(self,r=None,par=None,dnm=None):
        self.filt=str(r[6])
        if self.filt=='.': self.filt=np.nan
        if str(r[5])=='.': self.qual=np.nan
        else: self.qual=float(r[5])
        self.info_features(r[7])

    def output(self):
        o=[ self.n_alt,self.filt,self.qual,self.p_ar_max,self.p_ar_min,self.o_ar,
            self.p_dp_max,self.p_dp_min,self.o_dp,
            self.p_og_max,self.p_og_min,self.p_pg_max,self.p_pg_min,self.og,self.o_pg,
            self.p_gq_max,self.p_gq_min,self.o_gq]

        o+=self.info_feats
        return '\t'.join(map(str,o))

    def header(self):

        o = [   'nalt','filter','qual','parent_ar_max','parent_ar_min','offspring_ar',
                'parent_dp_max','parent_dp_min','offspring_dp',
                'parent_dnm_pl_max','parent_dnm_pl_min','parent_inh_pl_max','parent_inh_pl_min',
                'offspring_dnm_pl','offspring_inh_pl',
                'parent_gq_max','parent_gq_min','offspring_gq']
        o += self.info_keys
        
        return '\t'.join(map(str,o))