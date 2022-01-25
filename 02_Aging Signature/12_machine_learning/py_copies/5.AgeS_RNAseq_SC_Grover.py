#!/usr/bin/env python
# coding: utf-8

# Here is the sequence of operations to process single cell data and find relative ages of single cells
# Select only one data set in a time, skip the rest.
# Grover data have reverse order of appearance of old and young cells, also reversed aging score, which all together feels like an annotation flop. When adjusted to other sets in the same order (first young , then old cells), then it behaves equally to others. Make proper adjustments for this set. Now all script parameters are adjusted for three other data sets.
# 

# In[1]:


#import needed stuff
import math, os, pandas as pd
import numpy as np, matplotlib.pyplot as plt
#from sklearn import preprocessing as sp
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
items=getitems(f_in,"\t") #these are lists per row
cells=items[0] #these are samples names
#find samples and indexes for young and old
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


# In[9]:

#check how data look like
from scipy import sparse
all_n, all_L=[],[] #data as is or log-transformed
for line in items[1:]:
    values_n=[int(x) for x in line[1:]]
    values_L=[math.log2(int(x)+1) for x in line[1:]]
    all_n+=values_n
    all_L+=values_L
print("All data points", len(all_n))
sparse_matrix=sparse.find(all_n)
print("All non-zero", len(sparse_matrix[2]))
print("All zero", len(all_n)-len(sparse_matrix[2]))
print("Sparsity", (len(all_n)-len(sparse_matrix[2]))/len(all_n))
#%% visualise the data if interested
plt.figure(figsize=[15,5]) #adjust the figure height and width to you browser screen
plt.subplot(121)
plt.title("Original data")
plt.hist(all_n)
plt.subplot(122)
plt.hist(all_L)
plt.title("Log2 transformed data")


# If the figure on the left looks like one lonely block at zero, it means data are not log-transformed.
# If you see a shoulder in some range next to the block, it is likely log-transformed.

# In[10]:

#log-transform data if they are not
logged, genes=[],[]
for line in items[1:]:
    logged.append([math.log2(int(x)+1) for x in line[1:]])
    genes.append(line[0])
# In[11]:

#check how data are normalised
#data=np.array([x for x in items[1:]]).T
data=np.array(logged).T
totals, present=[],[]
for idx, d in enumerate(data):
    values=[float(x) for x in d]
 #   print("totals in sample ", idx, sum(values))
    totals.append(sum(values))
    present.append(len([x for x in values if x>0]))
plt.figure(figsize=[15,5]) 
plt.scatter(totals, present)
#Adjust title to the data you use!
plt.title("Log2 transformed Grover")
plt.xlabel("Total counts")
plt.ylabel("Genes present")
print("median totals", np.median(totals))


# In[12]:

#normalize if not yet done
new=[] #new will be normalized data
new.append(genes) #append gene names
t,p=[],[] #new totals and present
for i in range(len(data)):
    #here we normalize to the same median per totals per gene numbers
    #note a coefficient (4,5,6 or close)
    #you can leave it as it is, or change slightly to adjust normalized data to the same median as above
    values=[float(x)*5*present[i]/totals[i] for x in data[i]]
    new.append(values)
    t.append(sum(values))
    p.append(len([x for x in values if x>0]))
plt.figure(figsize=[15,5]) 
plt.scatter(t, p)
plt.title("Normalized")
plt.xlabel("Total counts")
plt.ylabel("Genes present")
print("median totals", np.median(t))


# In[13]: optional


#before you reduce the table to AS genes you can inquire for any gene in the set
nm="Flt3" #gene name of your choice
print(nm)
my_gene=[x for x in np.array(new).T if x[0]==nm][0][1:]
values=[float(x) for x in my_gene]
print("young cells", np.average(values[:max(youngID)]))
print("old cells", np.average(values[max(youngID):])) 


# In[14]:


#reduce the table only to the genes present in the AS list

selected=[x for x in np.array(new).T if x[0] in ASgenes]
data=[]
for line in selected:
    values=[float(x) for x in line[1:]]
    data.append(values)
genes=[x[0] for subset in selected for x in subset]
print("length selected", len(selected))
genes=[x[0] for x in selected]
print("genes", len(genes))
#print("selected", selected[0])
#%% optional check
if "Clca3a1" in genes:
    print("Clca3a1 still there")
else:
    print("Clca3a1 not in this set")
#%% 
# save data for clustering and predictions
df=pd.DataFrame(data=data, index=genes, columns=cells)
out_d="/Users/leonid/Documents/Publications_Presentations/2020/AgingSign/data/sources/selected_AS/"
df.to_csv(out_d+f_in[:-4]+"_selected.tsv", sep="\t", index=True)
# In[15]:


#draw heatmap
values=[]
for series in selected:
    v=[float(x) for x in series[1:]]
    values.append(v)
plt.figure(figsize=[15,5])
plt.imshow(np.array(values).T); #was values
plt.xticks(range(len(genes)), labels=genes, rotation='vertical');


# In[16]:


#plot normalized data and 
#check per gene (skip for RNAseq data)
import scipy.stats as ss
present, averages, totals, sigs=[],[],[],[]
for series in values:
    pre=[x for x in series if float(x)>0]
    #add t-test for young and old
    yo=[float(x) for x in series[:max(youngID)]]
    ol=[float(x) for x in series[max(youngID):]]
    stat, pval =ss.ttest_ind(yo, ol)
    print(values.index(series), genes[values.index(series)], stat, pval)
    if pval<0.05:
        sigs.append(genes[values.index(series)])
    ave=np.mean(series)
    tot=sum(series)
    present.append(len(pre))
    averages.append(ave)
    totals.append(tot)
for i in range(len(genes)):
    #can test fold diff and p-val for t-test
    #present as volcano, dot size by present
 #   print(genes[i],present[i],averages[i],totals[i])
    if genes[i] in sigs:
        plt.text(present[i],averages[i], genes[i])
plt.scatter(present, averages)
plt.xlabel("Present")
plt.ylabel("Averages")
#plt.hist(present)


# In[17]:


#make list of scores
new=np.array([x[1:] for x in selected])#.T
scores=[0]*len(new[0]) #generate empty list to collect scores per cell
#transpose back to genes in rows and normalize again per gene
#new=sp.normalize(data, axis=1, norm="l2")
#or not
print("check lengthes of genes and data", len(genes), len(new))

#for each column add values for pos and subtract for neg
for i in range(len(genes)):
    if genes[i] in pos:
        for j in range(len(scores)):
            scores[j]+=float(new[i][j])
    if genes[i] in neg:
        for j in range(len(scores)):
            scores[j]-=float(new[i][j])
print(scores[:10],"...\n", len(scores), min(scores), max(scores))
# In[19]:


#scale scores to 0-1 interval
sco=[(x-min(scores))/max(scores) for x in scores]
scores=sco
print("new scores", scores[:10],"...\n", min(scores), max(scores))


# In[20]:


#this is PCA for all data
X = new.T
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
print(len(new),len(genes))
for i in range(len(genes)):
    if genes[i]==gene_to_color:
        val=[float(x) for x in new[i]]
        break
val_to_color=[x/max(val) for x in val]
#you can check all data together
#print(len(pc1))
#print(len(pc2))
#print(len(val_to_color))


# In[24]:
#here we show PCA plots colored by AS or gene expression

thresh=max(youngID)
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
X=np.array(selected).T[1:] #was data
genes=np.array(selected).T[0] 
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


clf.predict(X) #show prediction per each cell


# In[31]:


plt.scatter(clf.predict(X), scores);


# In[ ]:


#I will make another script to do more of the predictions and classifications

