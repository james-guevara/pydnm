#!/usr/bin/env python3
import pandas as pd
from sklearn.externals import joblib
from pyDNM.Backend import get_path
import os,sys

def classify_dataframe(df, clf, ofh):
    df = df.dropna(axis=0,subset=df.columns[12:36])
    if df.empty: 
        print("Empty dataframe.")
        return 0
    X = df[df.columns[12:36]].values
    df["pred"] = clf.predict(X)
    df["prob"] = clf.predict_proba(X)[:,1]
    with open(ofh, "a") as f:
        df.to_csv(f, header=False)

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
    
    '''
    # Make separate dataframes for each type of variant: autosomal snp, autosomal indel,
    # X snp, Y snp, X indel, Y indel
    dfY = df.loc[df["chrom"] == "chrY"]
    dfX = df.loc[df["chrom"] == "chrX"]
    # Autosomal 
    df_auto= df.loc[(df["chrom"] != "chrY") & (df['chrom'] != 'chrX')]
    snv = df_auto.loc[(df_auto['ref'].str.len()==1) & (df_auto['alt'].str.len()==1)]
    indel = df_auto.loc[(df_auto['ref'].str.len()!=1) | (df_auto['alt'].str.len()!=1)]
    snv = snv.dropna(axis=0,subset=snv.columns[12:36])
    X_test = snv[snv.columns[12:36]].values
    indel = indel.dropna(axis=0,subset=indel.columns[12:36])
    X_test_indel = indel[indel.columns[12:36]].values
    if snv.empty:print("No SNVs to Classify")
    else:
            snv['pred'] = clf.predict(X_test)
            snv['prob'] = clf.predict_proba(X_test)[:, 1]
    if indel.empty:print("No INDELs to Classify")
    else:        
            indel['pred'] = clf_indels.predict(X_test_indel)
            indel['prob'] = clf_indels.predict_proba(X_test_indel)[:, 1]
    df_auto = pd.concat([snv,indel])
    
    if keep_fp == False:
            df_auto = df_auto.loc[df_auto['pred']==1]
    df_auto.to_csv(ofh,sep="\t",index=False)

    # sex chromosome variants dataframes (X chromosome)
    snvX = dfX.loc[(dfX['ref'].str.len()==1) & (dfX['alt'].str.len()==1)]
    indelX = dfX.loc[(dfX['ref'].str.len()!=1) | (dfX['alt'].str.len()!=1)]
    snvX = snvX.dropna(axis=0,subset=snvX.columns[12:36])
    X_test = snvX[snvX.columns[12:36]].values
    indelX = indelX.dropna(axis=0,subset=indelX.columns[12:36])
    X_test_indel = indelX[indelX.columns[12:36]].values
    if snvX.empty:print("No chrX SNVs to Classify")
    else:        
            snvX['pred'] = clf_chrX_snps.predict(X_test)
            snvX['prob'] = clf_chrX_snps.predict_proba(X_test)[:, 1]
    if indelX.empty:print("No chrY INDELs to Classify")
    else:        
            indelX['pred'] = clf_chrX_chrY_indels.predict(X_test_indel)
            indelX['prob'] = clf_chrX_chrY_indels.predict_proba(X_test_indel)[:, 1]
    dfX = pd.concat([snvX,indelX])
    if keep_fp == False:
            dfX = dfX.loc[dfX['pred']==1]
    dfX.to_csv(ofh,sep="\t",index=False,)
            
    # sex chromosome variants dataframes (Y chromosome)
    snvY = dfY.loc[(dfY['ref'].str.len()==1) & (dfY['alt'].str.len()==1)]
    indelY = dfY.loc[(dfY['ref'].str.len()!=1) | (dfY['alt'].str.len()!=1)]
    snvY = snvY.dropna(axis=0,subset=snv.columns[12:36])
    X_test = snvY[snvY.columns[12:36]].values
    indelY = indelY.dropna(axis=0,subset=snv.columns[12:36])
    X_test_indel = indelY[indelY.columns[12:36]].values
    if snvY.empty:print("No chrY SNVs to Classify")
    else:
            snvY['pred'] = clf_chrY_snps.predict(X_test)
            snvY['prob'] = clf_chrY_snps.predict_proba(X_test)[:, 1]
    if indelY.empty:print("No chrY INDELs to Classify")
    else:
            indelY['pred'] = clf_chrX_chrY_indels.predict(X_test_indel)
            indelY['prob'] = clf_chrX_chrY_indels.predict_proba(X_test_indel)[:, 1]
    #df = snv # TO DO: ADD IN THE INDELS!!!
    if not snvY.empty: dfY = pd.concat([snvY,indelY])
        if keep_fp == False:
                dfY = dfY.loc[dfY['pred']==1]
        dfY.to_csv(ofh,sep="\t",index=False)

    '''
