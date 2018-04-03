# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 15:21:55 2017

@author: Gayatri
"""

import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn.externals import joblib
def MODEL_TESTING(coverage_file,base_file,gc_file,mismatch_file,mfe_file,output_file):
    window=40
    w=window-1
    print('STARTING FEATURE EXTRACTION')
    with open(coverage_file) as f:
        lines = f.readlines()
    pos=[]
    coverage=[]
#    
    for items in lines:
        a=items.split()
        pos.append(int(a[1]))
        coverage.append(a[2])
    with open(base_file) as f:
        lines = f.readlines()
    nucleotide=[]
    for items in lines:
        a=items.split()
        nucleotide.append(a[2])
####################################################################
    loc10=[]
    gc_raw=[]
    char_raw=[]
    with open(gc_file) as f:
        lines = f.readlines()
    for items in lines:
        a=a=items.split()
        if a[3][0]!='N':
            loc10.append(a[1]+'-'+a[2])
            gc_raw.append(float(a[3][0:-1]))
            char_raw.append(a[0])
    loc11=[]
    mfe_raw=[]
    with open(mfe_file) as f:
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
    with open(mismatch_file) as f:
        lines = f.readlines()
    for items in lines:
        a=a=items.split()
        
        mismatch_whole.append(a[2])
####################################################################
    pos=list(map(int,pos))
    coverage=list(map(float,coverage))
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
#############################################################3
    for i in range(len(gc)):
        c20[i]=np.append(c20[i],gc[i])
        c20[i]=np.append(c20[i],mismatch[i])
        c20[i]=np.append(c20[i],mfe[i])
    X=c20
    
    #UNCOMMENT THIS PART IF YOU WANT OT OBSERVE YOUR FEATURES
#    df = pd.DataFrame(np.array(X))
#    df.to_csv("TrainingData.csv")
##########################################################################3
 
    LogReg  = joblib.load('LogReg_model.sav')
    poly = PolynomialFeatures(degree=2) # you can change this to suite your data
    Xa = poly.fit_transform(X)
    yhat_LogReg=LogReg.predict(Xa)
#    print 'AAAAAAAAAAAAAAAAAAAAAaaa'
#    print yhat_LogReg[0:2]
#    
    prob_LogReg= LogReg.predict_proba(Xa)
#    print prob_LogReg[0:2]
    p0=[]
    yhat_LogReg=[]
    for items in prob_LogReg:
        if items[0]>items[1]:
            yhat_LogReg.append(0)
            p0.append(items[0])
        else:
            yhat_LogReg.append(1)
            p0.append(items[1])
    
    SVM=joblib.load('SVM_model.sav')
    output=[]
#    yhat_SVM=SVM.predict(X)
    prob_SVM=SVM.predict_proba(X)
#    print SVM.predict(X)
    p=[]
    yhat_SVM=[]
    for items in prob_SVM:
        if items[0]>items[1]:
            yhat_SVM.append(0)
            p.append(items[0])
        else:
            yhat_SVM.append(1)
            p.append(items[1])
#    print prob_SVM[0:3]
    c=0
    for i in range(len(yhat_SVM)):
        
        if yhat_LogReg[i]==yhat_SVM[i]==1:
            c=c+1
            output.append( char[i] + '\t' + str(s1[i]) + '\t' +str(s2[i])+ '\t' + "%.2f" % ( (p[i]+p0[i]) /2.0) )
    print('Number of positives detected',len(output))
#    print output
    with open(output_file, 'w') as thefile:
        for item in output:
            thefile.write("%s\n" % item)
    thefile.close
    return 
#########################################################################
#coverage_file='positives_temp\\coverage.txt'       
#base_file='positives_temp\\refBase.txt' 
#gc_file='positives_temp\\gcContent.txt' 
#mismatch_file='positives_temp\\mismatch.txt'
#mfe_file='positives_temp\\mfe.txt'
#
#ans=MODEL_TESTING(coverage_file,base_file,gc_file,mismatch_file,mfe_file)
#coverage_file='negatives_temp\\coverage.txt'       
#base_file='negatives_temp\\refBase.txt' 
#gc_file='negatives_temp\\gcContent.txt' 
#mismatch_file='negatives_temp\\mismatch.txt'
#mfe_file='negatives_temp\\mfe.txt'
#
#ans=MODEL_TESTING(coverage_file,base_file,gc_file,mismatch_file,mfe_file)
##print(ans)
