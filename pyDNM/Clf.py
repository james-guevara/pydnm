#!/usr/bin/env python3
import pandas as pd
from sklearn.externals import joblib
from pyDNM.Backend import get_path
import os,sys
def classify(ofh=None,keep_fp=None):
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
	df = pd.read_csv(ofh,sep="\t")
	
	dfY = df.loc[df["chrom"] == "chrY"]
	dfX = df.loc[df["chrom"] == "chrX"]
	# df_auto = df.loc[~df["chrom"].str.contains("X")]
	# df_auto = df_auto.loc[~df_auto["chrom"].str.contains("Y")]
	df_auto= df.loc[(df["chrom"] != "chrY") & (df['chrom'] != 'chrX')]

	#print(dfX)



	snv = df_auto.loc[(df_auto['ref'].str.len()==1) & (df_auto['alt'].str.len()==1)]
	indel = df_auto.loc[(df_auto['ref'].str.len()!=1) | (df_auto['alt'].str.len()!=1)]
	#print(indel)
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
        #df = snv # TO DO: ADD IN THE INDELS!!!
	df_auto = pd.concat([snv,indel])
	#print(df_auto)
	#snv = snv.loc[snv["pred"]==1]
	#snv.to_csv(ofh,sep="\t",index=False)
	#snv["pos_0"] = snv["pos"]-1
	#snv_regions = snv[["chrom","pos_0","pos"]]
	#snv_regions.to_csv("tmp.bed",sep="\t",index=False,header=False)
	# df_auto = pd.concat([snv])
	if keep_fp == False:
		df_auto = df_auto.loc[df_auto['pred']==1]
	df_auto.to_csv(ofh,sep="\t",index=False)
	
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
	#df = snv # TO DO: ADD IN THE INDELS!!!
	dfX = pd.concat([snvX,indelX])
	if keep_fp == False:
		dfX = dfX.loc[dfX['pred']==1]
	dfX.to_csv(ofh,sep="\t",index=False,)
		
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
	dfY = pd.concat([snvY,indelY])
	if keep_fp == False:
		dfY = dfY.loc[dfY['pred']==1]
	dfY.to_csv(ofh,sep="\t",index=False)

