#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 15:30:43 2021
play with AS SC data, try to check what and how predict
@author: leonid
"""
#%%
import pandas as pd, matplotlib.pyplot as plt, numpy as np, os, seaborn as sns
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.metrics import confusion_matrix #, pairwise_distances
from sklearn.preprocessing import normalize
from sklearn.manifold import MDS
from sklearn.cluster import KMeans
#%%
os.chdir("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/12_machine_learning/selected_AS")
os.getcwd()
#%% load files with AS genes only, they are different length
files=[x for x in os.listdir() if "norm_log"in x]
fn=[]
for name in files:
    df=pd.read_csv(name, sep="\t", index_col=0)
    fn.append(df)
#%% find ages of cells from column names
ages=[]
for i in range(len(fn)):
    age=[]
    if i!=1:
        for name in fn[i].columns:
            if "young" in name:
                age.append(0)
            else:
                age.append(1)
    else:
        for name in fn[i].columns:
            if "4m" in name:
                age.append(0)
            else:
                age.append(1)
    ages.append(age)
#%% scale per gene
scaled=[]
for data in fn:
    new=normalize(data, norm='l2', axis=1) #for fows axis=1
    df=pd.DataFrame(new, columns=data.columns, index=data.index)
    scaled.append(df)
#%% make a heatmap if needed
# sns.heatmap(scaled[2],xticklabels=False, yticklabels=False)
# plt.title("File2")
# plt.ylabel("Genes")
# plt.xlabel("Cells")
#%% collect best predictors either via ADA or Forest
#apply scaled data
fn=scaled
all_predictors, ada_sorts=[],[]
for i in range(len(fn)):
    sorts=[]
    X=fn[i].T #was data
    y=ages[i]
    genes=fn[i].index
    #use one of two possible classifiers
    #clf = RandomForestClassifier(max_depth=2, random_state=0)
    clf=AdaBoostClassifier(random_state=0)
    clf.fit(X, y)
    print("File", files[i][:10], "score:", clf.score(X, y))
    imp=clf.feature_importances_
    sort=sorted(zip(imp, genes), reverse=True)
  #  print("Random forest age predictors")
    for i in range(len(sort)):
        if sort[i][0]>0:
     #       print(i, sort[i])
            all_predictors.append(sort[i][1]) 
            sorts.append(sort[i])
        else:
            break
    ada_sorts.append(sorts)
    y_predict= clf.predict(X) #show prediction per each cell
    print(confusion_matrix(y, y_predict)) #y_true, y_pred)
#%% make table of the best predictors from ada_sorts
dfs=[]
for data in ada_sorts:
    freqs=[x[0] for x in data]
    genes=[x[1] for x in data]
    df=pd.DataFrame(freqs, index=genes )
  #  df.columns=["freq"]
  #  print(df)
    dfs.append(df)
ada_best=pd.concat(dfs, axis=1, join="outer")
ada_best.columns=["Ma","Ki","Ko","Gr"]
ada_best.fillna(0, inplace=True)
ada_best["sum"]=ada_best["Ma"]+ada_best["Ki"]+ada_best["Ko"]+ada_best["Gr"]

srt=round(ada_best.sort_values(['sum'], ascending=False),3)
print(srt.head(n=50))
#%% save this table (one for all ADA or Forest scans)
# srt.to_csv("Best_RNDforest.tsv", sep="\t")
#%% check occurrence of age predictors in sets
#this part is optional
Age_predictors=[]
predictors=set(all_predictors)
counts=[]
for name in predictors:
    counts.append(all_predictors.count(name))
sort=sorted(zip(counts, predictors), reverse=True)
for i in range(len(sort)):
    print(i, sort[i][0], sort[i][1])
    if sort[i][0]>3:
        Age_predictors.append(sort[i][1])
    if sort[i][0]==1:
        break
print("Age predictors", Age_predictors)
#%% reduce the table only to the genes present in the reduced AS list or other list
ADA_list= ['Nupr1','Sult1a1', 'Fyb']
RND_list= ['Tmem176a', 'Tgm2', 'Sult1a1', 'Selp', 'Ramp2', 'Prtn3', 'Prcp', 
'Ppp1r16b', 'Nupr1', 'Mef2c', 'Mcm5', 'Ly6e', 'Lpl', 'Jun', 'Exoc6b', 'Dnmt3b', 
'Dhrs3', 'Clec1a', 'Aldh1a1']
best20=set(ADA_list+RND_list)
best5=ADA_list+['Selp', 'Clec1a'] #best ADA + next best from RND
best2=['Nupr1','Selp']
best1=['Selp']

red=[]
for data in fn:
    selected=[x for x in data.index if x in best5] #Age_predictors]
    reduced=data.loc[selected]
    red.append(reduced)
#%%
#check genes
#all_genes=[]
#for i in range(len(red)):
#    for item in red[i].index:
#        all_genes.append(item)
#for i, item in enumerate(set(all_genes)):
#    print(i, item, all_genes.count(item))
#%% Mt1 is missing on one set
#Age_predictors.remove('Mt1')
#%%
# now I can use one set for train, other for predict
names, result=[],[]
for ID in range(len(red)):
    rs=[]
    print("Train_by_set", ID, files[ID][:8])
    X_train, y_train=red[ID].T,ages[ID]
    for i in range(len(ages)):   
        X_pred, y_pred=red[i].T,ages[i]
        #from sklearn.datasets import make_classification
        #clf = RandomForestClassifier(max_depth=2, random_state=0)
        clf=AdaBoostClassifier(random_state=0)
        clf.fit(X_train, y_train)
        y_predict= clf.predict(X_pred) #show prediction per each cell
        mtx=confusion_matrix(ages[i], y_predict) #y_true, y_pred)
      #  right=mtx[0][0]+mtx[1][1]
      #  wrong=mtx[0][1]+mtx[1][0]
      #  success=right/(right+wrong)
        print(round(clf.score(X_pred,y_pred),3))#, round(success,3))
        rs.append(round(clf.score(X_pred,y_pred),4))
    names.append(files[ID][:8])
    result.append(rs)
#%% print resulting table
print("best20,  RNDforest")
res_df=pd.DataFrame(result, columns=names)
print(res_df)
res_df.to_csv("202101_outputs/Best5_scores_ADAboost.tsv", sep="\t")
# %% visualise single gene if selection is only one gene
plt.style.use('ggplot')
print(plt.style.available)
for i in range (len(red)):
    print(len(red[i].values[0]),len(ages[i]))
    plt.subplot(2,2,i+1)
    plt.title(files[i][:8])
    ID=i
    df0=pd.DataFrame(data={'expr': red[ID].values[0], 'age': ages[ID]})
    sns.stripplot(data=df0, y='expr', x='age')
# %% if more than one gene
lines=["X\tY\tages\tKMeans\tmismatches"]
plt.style.use('ggplot')

red=scaled #if you want the full set

for i in range(len(red)):
    ID, colorz=i, []
    #plt.title("MDS correlation by age")
    #matrix=pairwise_distances(red[ID].T, metric="correlation")
    embed = MDS(n_components=2)#, dissimilarity="precomputed")
    X_tr = embed.fit_transform(red[ID].values.T).T
    #K-means
    clusters = KMeans(n_clusters=2, random_state=0).fit_predict(X_tr.T)
    
    #check order of nuls and ones
    nuls, ones=[],[]
    for r in range(len(clusters)):
        if clusters[r]==0:
            nuls.append(r)
        else:
            ones.append(r)
    if sum(nuls)>sum(ones):
        new=[]
        for value in clusters:
            if value==0:
                new.append(1)
            else:
                new.append(0)
        clusters=new
    #colorz
    
    colorz,right,wrong=[],0,0
    for j in range(len(ages[ID])):
        if ages[ID][j]!=clusters[j]:
            colorz.append(2)
            wrong+=1
        else:
            colorz.append(ages[ID][j])
            right+=1
    success=right/(right+wrong)
    #print out data for dots
    print("File", files[i][:8], round(success,3))
    lines.append("File"+files[i]+str(round(success,3)))
    for ii in range(len(X_tr[0])):
        line=str(X_tr[0][ii])+"\t"+str(X_tr[1][ii])+"\t"+ str(ages[i][ii])+"\t"+ str(clusters[ii])+"\t"+\
            str(colorz[ii])
        lines.append(line)
      #  print(X_tr[0][ii],X_tr[1][ii],ages[i][ii],clusters[ii],colorz[ii])
#if using full set (red=scaled) then use ages[i] for coloring  
    plt.subplot(2,2,i+1)
    plt.title(files[i][:8])
    plt.scatter(X_tr[0],X_tr[1], 
                c=colorz, #for small subsets
            #    c=ages[i], #for full set
                edgecolor="black")
#%%
f_out=open("202101_outputs/Allgenes_MDS_ADAboost.tsv","w")
f_out.writelines("Best20 genes, RNDforest scaled\n")
for line in lines:
    f_out.writelines(line+"\n")
f_out.close()
# %% 
# clusters = KMeans(n_clusters=2, random_state=0).fit_predict(X_tr.T)
# plt.scatter(X_tr[0],X_tr[1], c=clusters, edgecolor="black")
# mtx=confusion_matrix(ages[ID], clusters)
# right=mtx[0][0]+mtx[1][1]
# wrong=mtx[0][1]+mtx[1][0]
# success=max(right,wrong)/(right+wrong)
# print(success)
# %%
#add color for mismatching assignments
# colorz=[]
# for i in range(len(ages[ID])):
#     if ages[ID][i]!=clusters[i]:
#         colorz.append(3)
#     else:
#         colorz.append(ages[ID][i])
# plt.scatter(X_tr[0],X_tr[1], c=colorz, edgecolor="black")
# %%
