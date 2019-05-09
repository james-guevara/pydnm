#!/usr/bin/env python3 
import pandas as pd
from sklearn.externals import joblib
from pyDNM.Backend import get_path
import os,sys

def classify_dataframe(df, clf, ofh,pyDNM_header=False, mode="a"):
    df = df.dropna(axis=0,subset=df.columns[12:36])
    if df.empty: 
        print("Empty dataframe.")
        return 0
    X = df[df.columns[12:36]].values
    df["pred"] = clf.predict(X)
    df["prob"] = clf.predict_proba(X)[:,1]
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

def classify(ofh=None,keep_fp=None,pseud=None):
    # fam 
    fam_fh = "/home/j3guevar/pydnm_src/pydnm/reach_ssc1-4.fam"
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
    df.to_csv("bleh.txt",sep="\t",header=True,index=False)
    sys.exit()

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

     
    classify_dataframe(df_autosomal_SNV,clf,ofh,True,"w")
    classify_dataframe(df_autosomal_indel,clf_indels,ofh)
    classify_dataframe(df_female_X_SNV,clf,ofh)
    classify_dataframe(df_female_X_indel,clf_indels,ofh)
    classify_dataframe(df_male_nonPAR_X_SNV,clf_chrX_snps,ofh)
    classify_dataframe(df_male_nonPAR_Y_SNV,clf_chrY_snps,ofh)    
    classify_dataframe(df_male_nonPAR_X_indel,clf_chrX_chrY_indels,ofh)    
    classify_dataframe(df_male_nonPAR_Y_indel,clf_chrX_chrY_indels,ofh)
    classify_dataframe(df_male_PAR_X_SNV,clf,ofh)
    classify_dataframe(df_male_PAR_Y_SNV,clf,ofh)
    classify_dataframe(df_male_PAR_X_indel,clf_indels,ofh)
    classify_dataframe(df_male_PAR_Y_indel,clf_indels,ofh)
