#!/usr/bin/env python3 
import pandas as pd
from sklearn.externals import joblib
from pyDNM.Backend import get_path
import os,sys
from pysam import VariantFile
import pybedtools


def classify_dataframe(df, clf, ofh,pyDNM_header=False, mode="a",keep_fp=False):
    df = df.dropna(axis=0,subset=df.columns[12:36])
    if df.empty: 
        # print("Empty dataframe.")
        return 0
    X = df[df.columns[12:36]].values
    df["pred"] = clf.predict(X)
    df["prob"] = clf.predict_proba(X)[:,1]
    if keep_fp == False:
        df = df.loc[df["pred"] == 1]
    with open(ofh, mode) as f:
        df.to_csv(f, sep="\t",header=pyDNM_header, index=False)

def get_sex(fam_fh):
    fam = open(fam_fh, "r")
    fam_dict = {}
    for line in fam:
        linesplit = line.rstrip().split("\t")
        iid = linesplit[1]
        sex = linesplit[4]
        fam_dict[iid] = sex 
    df = pd.Series(fam_dict).to_frame("sex")
    df["iid"] = df.index
    df.reset_index(inplace=True)
    df.drop(columns=["index"],inplace=True)
    return df

def classify(ofh=None,keep_fp=None,pseud=None,vcf=None,make_bed=True,make_vcf=True,fam_fh=None):
    ofh_new = ofh + ".preds"
    
    # fam 
    df_fam = get_sex(fam_fh)
    pseud_chrX = pseud["chrX"]
    pseud_chrX_interval_one = pseud_chrX[0]
    pseud_chrX_interval_two = pseud_chrX[1]
    pseud_chrY = pseud["chrY"]
    pseud_chrY_interval_one = pseud_chrY[0]
    pseud_chrY_interval_two = pseud_chrY[1]
    # Get classifiers
    snv_clf = get_path()+'/pydnm.snv.clf.joblib'
    indels_clf = get_path()+'/pydnm.indels.clf.joblib'
    snv_chrX_clf=get_path()+'/chrX_training_snps.joblib'
    snv_chrY_clf= get_path()+'/chrY_training_snps.joblib'
    indels_chrX_chrY_clf= get_path()+'/chrX_chrY_training_indels.joblib'
    if not os.path.isfile(snv_clf):
            sys.stderr.write('FATAL ERROR: {} CLASSIFIER NOT FOUND\n'.format(snv_clf))
            sys.exit(1)
            
    clf = joblib.load(snv_clf)
    clf_indels = joblib.load(indels_clf)
    clf_chrX_snps = joblib.load(snv_chrX_clf)
    clf_chrY_snps = joblib.load(snv_chrY_clf)
    clf_chrX_chrY_indels = joblib.load(indels_chrX_chrY_clf)
    # Make dataframe from input pydnm file
    df = pd.read_csv(ofh,sep="\t")
    
    # Filter original dataframe
    # df = df.loc[(df["offspring_gt"] != "0/0")] 


    df = pd.merge(df, df_fam, on="iid")

    df_autosomal_SNV = df.loc[(df["chrom"] != "chrY") & (df["chrom"] != "chrX") & (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)]
    df_autosomal_indel = df.loc[(df["chrom"] != "chrY") & (df["chrom"] != "chrX") & ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1))]
    df_female_X_SNV = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "2") & (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)]
    df_female_X_indel = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "2") & ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1))]
    
    df_male_nonPAR_X_SNV = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "1") & (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1) & ~(df["pos"].between(pseud_chrX_interval_one[0],pseud_chrX_interval_one[1]) | df["pos"].between(pseud_chrX_interval_two[0], pseud_chrX_interval_two[1]))]
    df_male_nonPAR_Y_SNV = df.loc[(df["chrom"] == "chrY") & (df["sex"] == "1") & (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1) & ~(df["pos"].between(pseud_chrY_interval_one[0],pseud_chrY_interval_one[1]) | df["pos"].between(pseud_chrY_interval_two[0], pseud_chrY_interval_two[1]))]
    df_male_nonPAR_X_indel = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "1") & ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1)) & ~(df["pos"].between(pseud_chrX_interval_one[0],pseud_chrX_interval_one[1]) | df["pos"].between(pseud_chrX_interval_two[0], pseud_chrX_interval_two[1]))]
    df_male_nonPAR_Y_indel = df.loc[(df["chrom"] == "chrY") & (df["sex"] == "1") & ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1)) & ~(df["pos"].between(pseud_chrY_interval_one[0],pseud_chrY_interval_one[1]) | df["pos"].between(pseud_chrY_interval_two[0], pseud_chrY_interval_two[1]))]

    df_male_PAR_X_SNV = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "1") & (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1) & (df["pos"].between(pseud_chrX_interval_one[0],pseud_chrX_interval_one[1]) | df["pos"].between(pseud_chrX_interval_two[0], pseud_chrX_interval_two[1]))]
    df_male_PAR_Y_SNV = df.loc[(df["chrom"] == "chrY") & (df["sex"] == "1") & (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1) & (df["pos"].between(pseud_chrY_interval_one[0],pseud_chrY_interval_one[1]) | df["pos"].between(pseud_chrY_interval_two[0], pseud_chrY_interval_two[1]))]
    df_male_PAR_X_indel = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "1") & ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1)) & (df["pos"].between(pseud_chrX_interval_one[0],pseud_chrX_interval_one[1]) | df["pos"].between(pseud_chrX_interval_two[0], pseud_chrX_interval_two[1]))]
    df_male_PAR_Y_indel = df.loc[(df["chrom"] == "chrY") & (df["sex"] == "1") & ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1)) & (df["pos"].between(pseud_chrY_interval_one[0],pseud_chrY_interval_one[1]) | df["pos"].between(pseud_chrY_interval_two[0], pseud_chrY_interval_two[1]))]


    classify_dataframe(df_autosomal_SNV,clf,ofh_new,True,"w")
    classify_dataframe(df_autosomal_indel,clf_indels,ofh_new)
    classify_dataframe(df_female_X_SNV,clf,ofh_new)
    classify_dataframe(df_female_X_indel,clf_indels,ofh_new)
    classify_dataframe(df_male_nonPAR_X_SNV,clf_chrX_snps,ofh_new)
    classify_dataframe(df_male_nonPAR_Y_SNV,clf_chrY_snps,ofh_new)    
    classify_dataframe(df_male_nonPAR_X_indel,clf_chrX_chrY_indels,ofh_new)    
    classify_dataframe(df_male_nonPAR_Y_indel,clf_chrX_chrY_indels,ofh_new)
    classify_dataframe(df_male_PAR_X_SNV,clf,ofh_new)
    classify_dataframe(df_male_PAR_Y_SNV,clf,ofh_new)
    classify_dataframe(df_male_PAR_X_indel,clf_indels,ofh_new)
    classify_dataframe(df_male_PAR_Y_indel,clf_indels,ofh_new)
    
    if make_bed: 
        # ofb = make_output_bed(ofh_new)
        ofb = "TEST_SNV.bed"
        a = pybedtools.BedTool(ofb)
        b = pybedtools.BedTool(vcf)
        a_and_b = b.intersect(a, u=True, wa=True, header=False, output="a_and_b.pybed2.bed")

    # if make_vcf: make_output_vcf(vcf,ofh)

    


    
def make_output_bed(ofh):
    ofb = "TEST_SNV.bed"
    fout = open(ofb,"w")
    f = open(ofh,"r")
    f.readline()
    dnm_bed = []
    for line in f:
        linesplit = line.rstrip().split("\t")
        chrom = linesplit[0]
        pos = linesplit[1]
        ref = linesplit[3]
        alt = linesplit[4]
        if len(ref) != 1 or len(alt) != 1: continue
        iid = linesplit[5]
        pred = linesplit[-2]
        prob = linesplit[-1]
        pos_0 = str(int(pos)-1)
        pos_1 = pos
        ID_col = "{}:{}:{}:{}:{}:{}:{}".format(chrom,pos,ref,alt,iid,pred,prob)
        newline = "{}\t{}\t{}\t{}\n".format(chrom,pos_0,pos_1,ID_col)
        fout.write(newline)
        dnm_bed.append(newline)
    return ofb 



# def make_output_vcf(vcf,ofh):
#     ofv = "TEST.vcf"
#     vcf_in = VariantFile(vcf) 
#     print(vcf_in)
#     new_header = vcf_in.header
#     new_header.info.add("pyDNM_pred",1,"Float","pyDNM prediction")
#     new_header.info.add("pyDNM_prob",1,"Float","pyDNM probability")
#     # vcf_out = VariantFile("-", "w", header=vcf_in.header)
#     dnm_bed = set()
#     # for line in f:
#     #     linesplit = line.rstrip().split("\t")
#     #     chrom = linesplit[0]
#     #     pos = linesplit[1]
#     #     ref = linesplit[3]
#     #     alt = linesplit[4]
#     #     iid = linesplit[5]
#     #     pred = linesplit[-2]
#     #     prob = linesplit[-1]
#     #     if pred == "0": continue
#     #     newline = [chrom,pos,ref,alt]
#     #     dnm_bed.add(newline) 
#     for rec in vcf_in.fetch():
#         dnm = [rec.chrom,rec.pos,rec.ref,rec.alts[0]]
# 
