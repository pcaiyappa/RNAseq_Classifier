# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 15:21:55 2017

@author: Gayatri
"""
from sklearn.cross_validation import cross_val_score
from sklearn import svm
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import PolynomialFeatures
import pandas as pd 
from sklearn import metrics
from sklearn.externals import joblib
def MODEL_TRAINING(coverage_pos,base_pos,gc_pos,mismatch_pos,mfe_pos,coverage_neg,base_neg,gc_neg,mismatch_neg,mfe_neg):
    window=40
    w=window-1
    print('STARTING FEATURE EXTRACTION')
    with open(coverage_pos) as f:
        lines = f.readlines()
    pos=[]
    coverage=[]
#    
    for items in lines:
        a=items.split()
        pos.append(int(a[1]))
        coverage.append(a[2])
    with open(base_pos) as f:
        lines = f.readlines()
    nucleotide=[]
    for items in lines:
        a=items.split()
        nucleotide.append(a[2])
####################################################################
    loc10=[]
    gc_raw=[]
    char_raw=[]
    with open(gc_pos) as f:
        lines = f.readlines()
    for items in lines:
        a=a=items.split()
        if a[3][0]!='N':
            loc10.append(a[1]+'-'+a[2])
            gc_raw.append(float(a[3][0:-1]))
            char_raw.append(a[0])
    loc11=[]
    mfe_raw=[]
    with open(mfe_pos) as f:
        lines = f.readlines()
    for items in lines:
        a=a=items.split()
        if a[3][0]!='N':
            loc11.append(a[1]+'-'+a[2])
            mfe_raw.append(float(a[3]))
    loc1=[]
    gc1=[]
    mfe1=[]
    char1=[]
    for i in range(len(loc10)):
        if loc10[i] in loc11:
            index=loc11.index(loc10[i])
            loc1.append(loc10[i])
            gc1.append(gc_raw[i])
            mfe1.append(mfe_raw[i])
            char1.append(char_raw[i])
    mismatch_whole=[]
    with open(mismatch_pos) as f:
        lines = f.readlines()
    for items in lines:
        a=a=items.split()
        
        mismatch_whole.append(a[2])
####################################################################
    pos=list(map(int,pos))
    coverage=list(map(int,coverage))
    y=[]
    s1=[]
    s2=[]
    c20=[]
    loc2=[]
    gc=[]
    mfe=[]
    base=[]
    char=[]
    mismatch_list=[]
    for i in range(len(pos)):
        try:
            if int(pos[i+w])==int(pos[i])+w:
                s=str(pos[i]-1)+'-'+str(pos[i+w])
                if 'A' in nucleotide[i:i+window] or 'a' in nucleotide[i:i+window]:
                    if s in loc1:
                        index=loc1.index(s)
                        gc.append(gc1[index])
                        mfe.append(mfe1[index])
                        char.append(char1[index])
                        s1.append(pos[i])
                        s2.append(pos[i+w])
                        mismatch_list.append(mismatch_whole[i:i+window])
                        c20.append(np.array(coverage[i:i+window])/sum(coverage[i:i+window]))
                        base.append(nucleotide[i:i+window])
                        loc2.append(str(pos[i])+'-'+str(pos[i+w]))
                        y.append(1)

        except:
                z=0
    c=0
    mismatch=[]
    for m,b in zip(mismatch_list,base):
        temp=[]
        for x1,x2 in zip(m,b):
            if x2=='A' or x2=='a':
                temp.append(x1)
        if temp:   
            mismatch.append(max(temp))
        else:
            mismatch.append(0)

    print('Number of Positives',len(s1))
    p_nos=len(s1)

#########################################################################NEGATIVES
    with open(coverage_neg) as f:
        lines = f.readlines()

    pos=[]
    coverage=[]
    nucleotide=[]
    for items in lines:
        a=items.split()
        pos.append(int(a[1]))
        coverage.append(a[2])
    with open(base_neg) as f:
        lines = f.readlines()
    nucleotide=[]
    for items in lines:
        a=items.split()
        nucleotide.append(a[2])
##########################################################
    loc10=[]
    gc_raw=[]
    char_raw=[]
    with open(gc_neg) as f:
        lines = f.readlines()
    c=0
    for items in lines:
        a=a=items.split()
        c=c+1
        if a[3][0]!='N':
            loc10.append(a[1]+'-'+a[2])
            gc_raw.append(float(a[3][0:-1]))
            char_raw.append(a[0])
    loc11=[]
    mfe_raw=[]
    with open(mfe_neg) as f:
        lines = f.readlines()
    for items in lines:
        a=a=items.split()
        if a[3][0]!='N':
            loc11.append(a[1]+'-'+a[2])
            mfe_raw.append(float(a[3]))
    loc1=[]
    gc1=[]
    mfe1=[]
    char1=[]
    for i in range( len(loc10)):
        if loc10[i] in loc11:
            index=loc11.index(loc10[i])
            loc1.append(loc10[i])
            gc1.append(gc_raw[i])
            mfe1.append(mfe_raw[i])
            char1.append(char_raw[i])
    with open(mismatch_neg) as f:
        lines = f.readlines()
    for items in lines:
        a=a=items.split()
        
        mismatch_whole.append(a[2])
######################################################################
    pos=list(map(int,pos))
    coverage=list(map(int,coverage))
   

    for i in range(len(pos)):
        try:
            if int(pos[i+w])==int(pos[i])+w:
                s=str(pos[i]-1)+'-'+str(pos[i+w])
                if 'A' in nucleotide[i:i+window] or 'a' in nucleotide[i:i+window]:
                    if s in loc1:
                   
                        index=loc1.index(s)
                        gc.append(gc1[index])
                        mfe.append(mfe1[index])
                        char.append(char1[index])
                        s1.append(pos[i])
                        s2.append(pos[i+w])
                        mismatch_list.append(mismatch_whole[i:i+window])
                        c20.append(np.array(coverage[i:i+window])/sum(coverage[i:i+window]))
                        base.append(nucleotide[i:i+window])
                        loc2.append(str(pos[i])+'-'+str(pos[i+w]))
                        y.append(0)
        except:
            z=0
    for m,b in zip(mismatch_list,base):
        temp=[]
        for x1,x2 in zip(m,b):
            if x2=='A' or x2=='a':
                temp.append(x1)
        if temp:   
            mismatch.append(max(temp))
        else:
            mismatch.append(0)

    print('Number of Negatives',len(s1)-p_nos)

#############################################################################
    for i in range(len(gc)):
        c20[i]=np.append(c20[i],gc[i])
        c20[i]=np.append(c20[i],mismatch[i])
        c20[i]=np.append(c20[i],mfe[i])
    X=c20
    
    #UNCOMMENT THIS PART IF YOU WANT OT OBSERVE YOUR FEATURES
#    df = pd.DataFrame(np.array(X))
#    df.to_csv("TrainingData.csv")
    
#########################################################
    print( 'MODEL TRAINING LOGISTIC REGRESSION')
    
    
    poly = PolynomialFeatures(degree=2) # you can change this to suite your data
    Xa = poly.fit_transform(X)
    model = LogisticRegression(verbose=1,max_iter=500) # you can change the max_iters to suite your data
    model = model.fit(Xa, y)
    filename = 'LogReg_model.sav'
    joblib.dump(model, filename)
    
    # USING CROSS VALIDATION CHECK THE ACUURACY OF THE CLASSIFIER
   
    yhat=model.predict(Xa)
    print ( metrics.classification_report(y, yhat))
    print('\nMODEL TRAINING LOGISTIC REGRESSION accuracy')
    scores = cross_val_score(LogisticRegression(), Xa, y, scoring='accuracy', cv=10)
    print(scores.mean()*100, '%')
    print(yhat[0:10])
    print(yhat[200:210])
    print(yhat[1000:1010])
    print(yhat[1800:1810])
    print(yhat[1200:1210])
#    return
##################################################################################################################SVM
    print('\nMODEL TRAINING SVM accuracy')   
    degree_svm=3 # you can change this to suite your data
 
    clf = svm.SVC(kernel='poly',degree=degree_svm,probability=True, C = 1.0)
    clf= clf.fit(X, y)
    filename = 'SVM_model.sav'
    joblib.dump(clf, filename)
   
     #  CHECK THE ACUURACY OF THE CLASSIFIER
    clf.fit(X,y)
    yhat1=clf.predict_proba(X)
    #yhat1 gives the confidence intervals if you wish to see them
    yhat=clf.predict(X)
    print((1-sum(abs(y-yhat))/len(y))*100,'%')
    return
##############################################

coverage_pos='positives_temp\\coverage.txt'       
base_pos='positives_temp\\refBase.txt' 
gc_pos='positives_temp\\gcContent.txt' 
mismatch_pos='positives_temp\\mismatch.txt'
mfe_pos='positives_temp\\mfe.txt'
coverage_neg='negatives_temp\\coverage.txt'       
base_neg='negatives_temp\\refBase.txt' 
gc_neg='negatives_temp\\gcContent.txt' 
mismatch_neg='negatives_temp\\mismatch.txt'
mfe_neg='negatives_temp\\mfe.txt'
MODEL_TRAINING(coverage_pos,base_pos,gc_pos,mismatch_pos,mfe_pos,coverage_neg,base_neg,gc_neg,mismatch_neg,mfe_neg)
