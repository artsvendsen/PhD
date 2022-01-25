#!/usr/bin/env python
# coding: utf-8

# Here is the sequence of operations to process single cell data and find relative ages of single cells
# Select only one data set in a time, skip the rest.
# Grover data have reverse order of appearance of old and young cells, also reversed aging score, which all together feels like an annotation flop. When adjusted to other sets in the same order (first young , then old cells), then it behaves equally to others. Make proper adjustments for this set. Now all script parameters are adjusted for three other data sets.
# 

# In[1]:


#import needed stuff
import math, os, pandas as pd, seaborn as sns
import numpy as np, matplotlib.pyplot as plt, matplotlib as mpl
#from sklearn import preprocessing as sp
from scipy import sparse
from sklearn.decomposition import PCA
#I still like to use own function to get data from file
def getitems(fil, split):
    sfil=open(fil, "r"); lines=sfil.readlines(); sfil.close()
    items=[x.rstrip().replace("'","").replace('"','').split(split) for x  in lines]
    return(items)
print("Imports")


# In[3]:
"""
#First we load AL file and select AS genes. We have to keep separate negative and positive AS genes 
Positive will add to the score as they higher and negative will add as they are lower
"""
#read the file
dr = '/Users/leonid/Documents/Publications_Presentations/2020/AgingSign/data/aging_py/'
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
# In[4]:
# probably good place to set default directory
# define it where your source files are
os.chdir("/Users/leonid/Documents/Publications_Presentations/2020/AgingSign/data/sources/")
os.getcwd()
# In[5]:
f_in="Grover_count_filteredcorr.txt"
df_SC=pd.read_csv(f_in, sep="\t")
#%%
cells=df_SC.columns #these are samples names
#find samples and indexes for young and old, age is for color etc
youngID, oldID, age=[],[],[]
for i in range(len(cells)):
    if "young" in cells[i]:
        youngID.append(i)
        age.append(0)
    if "old" in cells[i]:
        oldID.append(i)
        age.append(1)
print("young", min(youngID),max(youngID)) #indexes for young cells
print("old", min(oldID),max(oldID)) #indexes for old cells
#%% redundant, age is for the same
#define groups
#groups=[]
#for name in df_SC.columns:
#    if "young" in name:
#        groups.append(0)
#    else:
#        groups.append(1)
#%% Check sparsity, a.k.a. how many zeros in the table
non_zero_points=sparse.find(df_SC)[2]
non_zero=len(non_zero_points) #len of all non-zero values found
all_points=len(df_SC)*len(df_SC.columns)
print("Sparcity", non_zero/all_points)

#%% visualise the data if interested
plt.figure(figsize=[15,5]) #adjust the figure height and width to you browser screen
plt.subplot(121)
plt.title("Original data")
plt.hist(non_zero_points)
plt.subplot(122)
L=[math.log2(x) for x in non_zero_points]
plt.hist(L)
plt.title("Log2 transformed data")
# If the figure on the left looks like one lonely block at zero, it means data are not log-transformed.
# If you see a shoulder in some range next to the block (like on the right side), it is likely log-transformed.
#%% define present and totals, show distribution
def get_pres_totl(df):
    pre,tot=[],[]
    for name in df.columns:
        total=sum(df[name])
        tot.append(total)
        present=np.count_nonzero(df[name])
        pre.append(present)
    return(pre,tot)
#type style.available to see possible styles
mpl.style.use('ggplot') #classic or default
pre,tot=get_pres_totl(df_SC)
plt.scatter(pre,tot, c=age, edgecolor="black")
#%% normalise and log transform in that order
medianT=np.median(tot)
medianP=np.median(pre)
for i in range(len(df_SC.columns)):
    data=df_SC[df_SC.columns[i]]
    step1=data*medianT/tot[i]
    step2=step1*pre[i]/medianP
    step3=np.log2(step2+1)
    df_SC[df_SC.columns[i]]=step3
pre,tot=get_pres_totl(df_SC)
plt.scatter(pre,tot, c=age, edgecolor="black")
# In[10]:
gene_col=df_SC.loc['Flt3']
yo,ol=[],[]
for i in range(len(age)):
    if age[i]==0:
        yo.append(gene_col[i])
    else:
        ol.append(gene_col[i])

print("young cells", np.average(yo))
print("old cells", np.average(ol) )


# In[14]:

#reduce the table only to the genes present in the AS list
selected=[x for x in df_SC.index if x in ASgenes]
data=df_SC.loc[selected]
#print("selected", selected[0])
#%% optional check
if "Clca3a1" in selected:
    print("Clca3a1 still there")
else:
    print("Clca3a1 not in this set")
#%% 
# save data for clustering and predictions

out_d="/Users/leonid/Documents/Publications_Presentations/2020/AgingSign/data/sources/selected_AS/"
data.to_csv(out_d+f_in[:-4]+"_sel_norm_log.tsv", sep="\t", index=True)
# In[15]:
#draw heatmap
sns.heatmap(data,square=True, yticklabels=True, xticklabels=False)

# In[16]:

#plot averages to present
#check per gene significance
import scipy.stats as ss
present, averages, totals, sigs=[],[],[],[]
for genes in data.index:
    series=data.loc[genes]
    pre=[x for x in series if float(x)>0]
    #add t-test for young and old
    y=[x for x in series.index if "young"in x]
    o=[x for x in series.index if "old"in x]
    yo=series[y]
    ol=series[o]
    stat, pval =ss.ttest_ind(yo, ol)
    print(genes, np.mean(yo),np.mean(ol), pval)
    if pval<0.000001:
        sigs.append(genes)
    ave=np.mean(series)
    tot=sum(series)
    present.append(len(pre))
    averages.append(ave)
    totals.append(tot)
for counter, gene in enumerate(data.index):
    #can test fold diff and p-val for t-test
    #present as volcano, dot size by present
 #   print(genes[i],present[i],averages[i],totals[i])
    if gene in sigs:
        plt.text(present[counter],averages[counter], gene)
plt.scatter(present, averages)
plt.xlabel("Present")
plt.ylabel("Averages")
#plt.hist(present)


# In[17]:
#make list of scores, lists of pos and neg are previously made
genes=data.index
scores=[0]*len(data.columns) #generate empty list to collect scores per cell

#for each column add values for pos and subtract for neg
for i in range(len(genes)):
    if genes[i] in pos:
        for j in range(len(scores)):
            scores[j]+=float(data.values[i][j])
    if genes[i] in neg:
        for j in range(len(scores)):
            scores[j]-=float(data.values[i][j])
print(scores[:10],"...\n", len(scores), min(scores), max(scores))
# In[19]:


#scale scores to 0-1 interval
sco=[(x-min(scores))/max(scores) for x in scores]
scores=sco
print("new scores", scores[:10],"...\n", min(scores), max(scores))


# In[20]:


#this is PCA for all data
X = data.T
pca = PCA(n_components=3)
pcomp=pca.fit_transform(X)
pc1=pcomp.T[0]
pc2=pcomp.T[1]
print(len(pc1), len(pc2))
print("Explained variance", pca.explained_variance_ratio_)


# In[21]:
#here we check the most contributing genes to component 1 and 2

pc1_genes=pca.components_[0]
print("len pc1", len(pc1_genes))
print("Genes in 1st component")
sort=sorted(zip(pc1_genes, genes), reverse=True)
for i in range(len(sort)):
    if abs(sort[i][0])>0.12:
        print(i, sort[i])
pc2_genes=pca.components_[1]
print("Genes in 2nd component")
sort=sorted(zip(pc2_genes, genes), reverse=True)
for i in range(len(sort)):
    if abs(sort[i][0])>0.12:
        print(i, sort[i])


# In[23]:


#try to pick some gene for coloring
gene_to_color= "Selp"
val=data.loc["Selp"]
val_to_color=[x/max(val) for x in val]
#you can check all data together
#print(len(pc1))
#print(len(pc2))
#print(len(val_to_color))


# In[24]:
#here we show PCA plots colored by AS or gene expression

thresh=min(oldID)
plt.figure(figsize=[15,5])
plt.subplot(121)
plt.title("AS scores")
plt.scatter(pc1[:thresh],pc2[:thresh],marker="v", edgecolor="black",
            cmap="Purples", 
            c=scores[:thresh])
plt.scatter(pc1[thresh:],pc2[thresh:],marker="o",edgecolor="black",
            cmap="Purples", 
            c=scores[thresh:])
plt.colorbar()
plt.clim(0,1)
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.subplot(122)
plt.scatter(pc1[:thresh],pc2[:thresh],marker="v",edgecolor="black",
            cmap="Purples", 
            c=val_to_color[:thresh])
plt.scatter(pc1[thresh:],pc2[thresh:],marker="o",edgecolor="black",
            cmap="Purples", 
            c=val_to_color[thresh:])
#for i in range(len(cells)):
 #   plt.text(pc1[i],pc2[i],cells[i][0:1])

plt.title(gene_to_color+ " scores")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.colorbar()
plt.clim(0,1)


# In[25]:
#more composite plot of the same

plt.figure(figsize=[15,10])
plt.subplot(221)
plt.title("Young")
#for i in range(0,48): 
 #   plt.text(pc1[i],pc2[i], cells[i][:3])
plt.scatter(pc1[0:thresh],pc2[0:thresh], c=scores[0:thresh])
plt.clim(0,1)
plt.colorbar()    
plt.xlim(-25,25)
plt.ylim(-15,15)
plt.subplot(222)
#for i in range(47,len(cells)):
 #   plt.text(pc1[i],pc2[i], cells[i][:3])
plt.scatter(pc1[thresh:],pc2[thresh:], c=scores[thresh:])
plt.clim(0,1)
plt.colorbar()
plt.xlim(-25,25)
plt.ylim(-15,15)
plt.title("Old")
plt.subplot(223)
plt.title("Mean AS in young and old")
plt.boxplot([scores[:thresh], scores[thresh:]])
plt.xticks([1,2], ["young", "old"])
plt.subplot(224)
for i in range(len(genes)):
    #can test fold diff and p-val for t-test
    #present as volcano, dot size by present
 #   print(genes[i],present[i],averages[i],totals[i])
    if genes[i] in sigs:
        plt.scatter(present[i],averages[i], c='red')
        #plt.text(present[i],averages[i], genes[i])
    else:
        plt.scatter(present[i], averages[i], c='gray')
plt.xlabel("Present")
plt.ylabel("Averages")
plt.title("Significant")
#plt.scatter(t, p)
plt.show()


# Note that aging score does not fully predict (or separate) young and old stem cells. There are some cells in young which are as old as some old cells (and other way around). For proper prediction we can use some machine learning algorythms, like for instance Random Forest.   

# In[26]:


"""extra if you want it. Find how successful 
is prediction based on any kind of selected data
first, define which data to use and give indexes of ages
"""
#run random forest
X=data.T #was data
y=age
from sklearn.ensemble import RandomForestClassifier
#from sklearn.datasets import make_classification
clf = RandomForestClassifier(max_depth=2, random_state=0)
clf.fit(X, y)
print("success of prediction, score:", clf.score(X, y))
imp=clf.feature_importances_
sort=sorted(zip(imp, genes), reverse=True)
print("Random forest age predictors")
for i in range(len(sort))[:10]:
    if sort[i][0]>0:
        print(i, sort[i])


# In[28]:


y_predict= clf.predict(X) #show prediction per each cell

from sklearn.metrics import confusion_matrix
print(confusion_matrix(y, y_predict)) #y_true, y_pred)
# In[31]:

sns.set_theme(style="darkgrid")
sns.stripplot(y=scores, x=(age), jitter=0.1)
#plt.scatter(clf.predict(X), scores);


# In[ ]:


#I will make another script to do more of the predictions and classifications

