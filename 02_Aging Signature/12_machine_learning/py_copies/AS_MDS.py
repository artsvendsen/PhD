#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 13:45:47 2021
show AS tables in MDS, check for their natural separation
@author: leonid
"""
from sklearn.manifold import MDS
import pandas as pd, matplotlib.pyplot as plt, numpy as np, os, seaborn as sns
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.metrics import confusion_matrix, pairwise_distances
from sklearn.preprocessing import normalize
from sklearn.cluster import KMeans
os.chdir("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/12_machine_learning/selected_AS")
os.getcwd()
#%%
files=[x for x in os.listdir() if "norm_log"in x]
fn=[]
for name in files:
    df=pd.read_csv(name, sep="\t", index_col=0)
    fn.append(df)
#%%
#load AS list
def getitems(fil, split):
    sfil=open(fil, "r"); lines=sfil.readlines(); sfil.close()
    items=[x.rstrip().replace("'","").replace('"','').split(split) for x  in lines]
    return(items)
#read the file
dr = '/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/09_single_cell/20200831_gdrive_AS_Plots/sources/'
fi= 'Aging_List_REAN_Feb6.txt'
items=getitems(dr+fi, "\t")[1:]
#get AS genes, those with citation >3
ASlines=[x for x in items if float(x[1])>3]
print("ASlines", ASlines[:10],"...")
#get list of negative genes, those which fold change is negative
neg=[x[0] for x in ASlines if float(x[2])<0]
print("neg", neg[:5],"...")
#similarly, get positive AS genes
pos=[x[0] for x in ASlines if float(x[2])>0]
print("pos", pos[:5],"...")
ASgenes=neg+pos
nm="Flt3"#Clca3a1
if nm in ASgenes:   
    print(nm, " in ASgenes")
#%% make ages for all datasets
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
#%%
#scale per gene
scaled=[]
for data in fn:
    new=normalize(data, norm='l2', axis=1) #for fows axis=1
    df=pd.DataFrame(new, columns=data.columns, index=data.index)
    scaled.append(df)
#%%
ID=0
data=fn[ID]
colorz=ages[ID]

#%%
plt.style.use('ggplot')
print(plt.style.available)
#%%
plt.subplot(221)
#plt.subplot(221)
plt.title("MDS correlation by age")
matrix=pairwise_distances(data.T, metric="correlation")
embed = MDS(n_components=2, dissimilarity="precomputed")
X_tr = embed.fit_transform(matrix).T
plt.scatter(X_tr[0],X_tr[1], c=ages[ID], edgecolor="black")
clusters = KMeans(n_clusters=2, random_state=0).fit_predict(data.T)
print("Ages", colorz)
print("Clusters", clusters)
plt.subplot(222)
plt.title("MDS by KMeans")
plt.scatter(X_tr[0],X_tr[1], c=clusters, edgecolor="black")
#Compute clustering and transform X to cluster-distance space
plt.subplot(223)
plt.title("K-means by ages")
kmeans_XY= KMeans(n_clusters=2, random_state=0).fit_transform(data.T)
plt.scatter(kmeans_XY.T[0],kmeans_XY.T[1], c=ages[ID], edgecolor="black")#ages[ID])
#options for clustering, color by ages
plt.subplot(224)
plt.title("K-means by KMeans")
#kmeans_XY= KMeans(n_clusters=2, random_state=0).fit_transform(data.T)
plt.scatter(kmeans_XY.T[0],kmeans_XY.T[1], c=clusters, edgecolor="black")
#%% #efficiency recalculated
print("Data file", files[ID])
print(confusion_matrix(ages[ID], clusters))
results=confusion_matrix(ages[ID],clusters)
correct=results[0][0]+results[1][1]
incorr=results[0][1]+results[1][0]
print("efficiency of prediction", max(correct,incorr)/(correct+incorr))


#%% try clustering only MDS coordinates
#mark mismatching cells
colorz=[]
for i in range(len(clusters)):
    if clusters[i]!=ages[ID][i]:
        colorz.append(3)
    else:
        colorz.append(ages[ID][i])
plt.title("MDS correlation by KMeans")
clusters = KMeans(n_clusters=2, random_state=0).fit_predict(X_tr.T)
plt.scatter(X_tr[0],X_tr[1], c=colorz, edgecolor="black", s=70);

# %%
