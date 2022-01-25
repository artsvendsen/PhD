"""
calculate aging index in expression data
Based on AS version 6Feb2020
by L.Bystrykh dec 2019
"""
#%%
#import needed stuff
import math
import numpy as np, matplotlib.pyplot as plt
#from sklearn import preprocessing as sp
from sklearn.decomposition import PCA
#I still like to use own function to get data from file
def getitems(fil, split):
    sfil=open(fil, "r"); lines=sfil.readlines(); sfil.close()
    items=[x.rstrip().replace("'","").replace('"','').split(split) for x  in lines]
    return(items)

#%%
"""
#First we load AL file and select AS genes. We have to keep separate negative and positive AS genes 
Positive will add to the score as they higher and negative will add as they are lower
"""
#read the file
dr = '/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/09_single_cell/20200831_gdrive_AS_Plots/sources/'
fi= 'Aging_List_REAN_Feb6.txt'
items=getitems(dr+fi, "\t")[1:]
#get AS genes, those with citation >3
ASlines=[x for x in items if float(x[1])>3]
print(ASlines)
#get list of negative genes, those which fold change is negative
neg=[x[0] for x in ASlines if float(x[2])<0]
print("neg", neg)
#similarly, get positive AS genes
pos=[x[0] for x in ASlines if float(x[2])>0]
print("pos", pos)
ASgenes=neg+pos
if "Clca3a1" in ASgenes:
    print("Clca3a1 in ASgenes")
print(len(ASgenes))
#%% #if you want data from Grover
dir="/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/09_single_cell/20200831_gdrive_AS_Plots/sources/"
f_in="Grover_count_filteredcorr.txt"
items=getitems(dir+f_in,"\t") #these are lists per row
cells=items[0][1:] #these are samples names

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
#%%
#try data Kowalczyk, they are in two files
f_in="/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/09_single_cell/20200831_gdrive_AS_Plots/sources/GSE59114_C57BL6_merged.csv"
file="GSE59114_C57BL6_merged.csv" 
items=getitems(f_in,"\t")
#clean from junk
cells=items[0][1:] #cell names
#check age groups
youngID, oldID, age=[],[],[]
for i in range(len(cells)):
    if "young" in cells[i]:
        youngID.append(i) #collect indexes for young cells
        age.append(0)
    if "old" in cells[i]:
        oldID.append(i) #collect indexes for old cells
        age.append(1)
print("young cells indexes", min(youngID),max(youngID))
print("old cells indexes", min(oldID),max(oldID))
#%%
#try Kirshner data, only 4m (young) and 18m (old). 
#Middle time point looks strange, better not to use
f_in="/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/09_single_cell/20200831_gdrive_AS_Plots/sources/Kirshner_2t_merged.csv"
items=getitems(f_in,"\t")
#clean from junk
cells=items[0][1:] #cell names
#check age groups
youngID, oldID, age=[],[],[]
for i in range(len(cells)):
    if "_4m" in cells[i]:
        youngID.append(i) #collect indexes for young cells
        age.append(0)
    if "_18m" in cells[i]:
        oldID.append(i) #collect indexes for old cells
        age.append(1)
print("young cells indexes", min(youngID),max(youngID))
print("old cells indexes", min(oldID),max(oldID))
#%% Mann SC data
#f_in="GSE100426_LT_filtered_Mann.txt"
dir="/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/09_single_cell/20200831_gdrive_AS_Plots/sources/"
f_in="GSE100426_SC_Mann_nostLTfilt.csv"
items=getitems(dir+f_in,"\t")
#clean from junk
cells=items[0][1:] #cell names
#check age groups
youngID, oldID, age=[],[],[]
for i in range(len(cells)):
    if "young" in cells[i]:
        youngID.append(i) #collect indexes for young cells
        age.append(0)
    if "old" in cells[i]:
        oldID.append(i) #collect indexes for old cells
        age.append(1)
print("young cells indexes", min(youngID),max(youngID))
print("old cells indexes", min(oldID),max(oldID))
#%%
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
# plt.subplot(121)
# plt.title("Original data")
# plt.hist(all_n)
# plt.subplot(122)
# plt.hist(all_L)
# plt.title("Log2 transformed data")
#%%
#log-transform data if they are not
logged, genes=[],[]
for line in items[1:]:
    logged.append([math.log2(int(x)+1) for x in line[1:]])
    genes.append(line[0])
    
#%%
#check how data are normalised
#data=np.array([x for x in items[1:]]).T
data=np.array(logged).T
totals, present=[],[]
for idx, d in enumerate(data):
    values=[float(x) for x in d]
 #   print("totals in sample ", idx, sum(values))
    totals.append(sum(values))
    present.append(len([x for x in values if x>0]))
# plt.scatter(totals[:min(oldID)], present[:min(oldID)])
# plt.scatter(totals[min(oldID):], present[min(oldID):])

# plt.title("Before normalization")
# plt.xlabel("Total counts")
# plt.ylabel("Genes present")
print("median totals", np.median(totals))
#%%
#normalize if not yet done
new=[] #new will be normalized data
new.append(genes) #append gene names
t,p=[],[] #new totals and present
for i in range(len(data)):
    #here we normalize to the same median per totals per gene numbers
    values=[float(x)*5*present[i]/totals[i] for x in data[i]]
    new.append(values)
    t.append(sum(values))
    p.append(len([x for x in values if x>0]))
# plt.scatter(t, p)
# plt.title("Normalized")
# plt.xlabel("Total counts")
# plt.ylabel("Genes present")

#%%
#reduce the table only to the genes present in the AS list
selected=[x for x in np.array(new).T if x[0] in ASgenes]
print("length selected", len(selected))
genes=[x[0] for x in selected]
print("genes", len(genes))
#print("selected", selected[0])
if "Clca3a1" in genes:
    print("Clca3a1 still there")
else:
    print("Clca3a1 not in this set")
#%%
#draw heatmap
values=[]
for series in selected:
    v=[float(x) for x in series[1:]]
    values.append(v)
plt.imshow(np.array(values).T) #was values
plt.xticks(range(len(genes)), labels=genes, rotation='vertical')

# Write out the values for this plot:
age = ["yng"] * len(youngID) + ["old"] * int(len(values)-len(youngID))
#age = ["old"] * int(len(values)-len(youngID)) + ["yng"] * thresh #for Mann dataset only
f_out=open("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/09_single_cell/normalized_Mann.txt","w")
f_out.writelines("genes\tvalue\n")
f_out.writelines(","+str(cells)+"\n")
for i in range(len(genes)):
    f_out.writelines(str(genes[i])+"\t"+str(values[i])+"\n")
f_out.close()

#%%
#plot normalized data and 
#check per gene (skip for RNAseq data)
import scipy.stats as ss
present, averages, totals, sigs, pval2=[],[],[],[],[]
print("x", "y", "gene", "pval")
f_out=open("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/09_single_cell/02_significance_Mann.txt","w")
f_out.writelines("x\ty\tgene\tpval\n")
for series in values:
    pre=[x for x in series if float(x)>0]
    #add t-test for young and old
    yo=[float(x) for x in series[:max(youngID)]]
    ol=[float(x) for x in series[max(youngID):]]
    stat, pval =ss.ttest_ind(yo, ol)
    # print(values.index(series), genes[values.index(series)], stat, pval)
    pval2.append(pval)
    if pval<0.001:
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
    #print(present[i],averages[i], genes[i], pval2[i])
    f_out.writelines(str(present[i])+"\t"+str(averages[i])+"\t"+genes[i]+"\t"+str(pval2[i])+"\n")
f_out.close()

plt.scatter(present, averages)
plt.xlabel("Present")
plt.ylabel("Averages")
#plt.hist(present)
#%%
#make empty list of scores
new=np.array([x[1:] for x in selected])#.T
scores=[0]*len(new[0]) #generate empty list to collect scores per cell

print("check lengthes of genes and data", len(genes), len(new))
#%%
#for each column add values for pos and subtract for neg
for i in range(len(genes)):
    if genes[i] in pos:
        for j in range(len(scores)):
            scores[j]+=float(new[i][j])
    if genes[i] in neg:
        for j in range(len(scores)):
            scores[j]-=float(new[i][j])
print(scores,"\n", len(scores), min(scores), max(scores))
#%%
#scale scores to 0-1 interval
sco=[(x-min(scores))/max(scores) for x in scores]
scores=sco
print("new scores", scores, "\n", min(scores), max(scores))
#%%
#this is PCA for all data
X = new.T #PCA works horizontally, cells should be in rows
pca = PCA(n_components=4)
pcomp=pca.fit_transform(X).T
print(len(pcomp),len(pcomp[0]))
pc1=pcomp[0]
pc2=pcomp[1]
print("PCA explained variance", pca.explained_variance_ratio_)
#%% this part is not working yet
#check the biggest factors in each component
#size1=pca.explained_variance_ratio_
#print("Size1", size1)
#size2=float(pca.explained_variance_ratio_[1])
#sizes=[size1, size2]
#print(len(genes))
#print(len(size1))
#zipped=zip(genes, size1)
#sort=sorted(zipped, reversed=True)
#for i in range(10):
#    print(sort)
#pc1=pca.components_[0] #plot coordinates
#pc2=pca.components_[1]
#%%
#try to pick some gene for coloring
gene_to_color="Selp"
print(len(new),len(genes))
for i in range(len(genes)):
    if genes[i]==gene_to_color:
        val=[float(x) for x in new[i]]
     #   break
val_to_color=[x/max(val) for x in val]
#%%
#you can check all data together
print(len(pc1))
print(len(pc2))
print(len(val_to_color))
#%%
#PCA AND AGING SIGNATURE SCORES
#thresh=max(youngID)
thresh=min(youngID) #for Mann dataset only
plt.subplot(121)
plt.title("AS scores")
plt.scatter(pc1[:thresh],pc2[:thresh],marker="v",
            cmap="Purples", 
            c=scores[:thresh])
plt.scatter(pc1[thresh:],pc2[thresh:],marker="o",
            cmap="Purples", 
            c=scores[thresh:])
plt.colorbar()
plt.clim(0,1)
plt.xlabel("PC1 ")
plt.ylabel("PC2 ")

# Write out the values for this plot:
#age = ["yng"] * thresh + ["old"] * int(len(pc1)-thresh)
age = ["old"] * int(len(pc1)-thresh) + ["yng"] * thresh #for Mann dataset only
f_out=open("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/09_single_cell/01_PCA_Mann_AS_scores.txt","w")
f_out.writelines("PC1\tPC2\tscore\tage\n")

for i in range(len(pc1)):
    #print(i)
    f_out.writelines(str(pc1[i])+"\t"+str(pc2[i])+"\t"+str(scores[i])+"\t"+str(age[i])+"\n")
f_out.close()


# plt.subplot(122)
# plt.scatter(pc1[:thresh],pc2[:thresh],marker="v",
#             cmap="Purples", 
#             c=val_to_color[:thresh])
# plt.scatter(pc1[thresh:],pc2[thresh:],marker="o",
#             cmap="Purples", 
#             c=val_to_color[thresh:])
# #for i in range(len(cells)):
#  #   plt.text(pc1[i],pc2[i],cells[i][0:1])
# plt.title("Selp scores")
# plt.xlabel("PC1 ")
# plt.ylabel("PC2 ")
# plt.colorbar()
# plt.clim(0,1)
#%%
# plt.title("Grover SC")
# #plt.colorbar()

# plt.subplot(221)
# plt.title("Young")
# #for i in range(0,48): 
#  #   plt.text(pc1[i],pc2[i], cells[i][:3])
# plt.scatter(pc1[0:thresh],pc2[0:thresh], c=scores[0:thresh])
# plt.clim(0,1)
# plt.colorbar()    
# plt.xlim(-25,25)
# plt.ylim(-25,25)
# plt.subplot(222)
# #for i in range(47,len(cells)):
#  #   plt.text(pc1[i],pc2[i], cells[i][:3])
# plt.scatter(pc1[thresh:],pc2[thresh:], c=scores[thresh:])
# plt.clim(0,1)
# plt.colorbar()
# plt.xlim(-25,25)
# plt.ylim(-25,25)
# plt.title("Old")
# plt.subplot(223)
# plt.title("Mean AS in young and old")
# plt.boxplot([scores[:thresh], scores[thresh:]])
# plt.xticks([1,2], ["young", "old"])
# plt.subplot(224)
# for i in range(len(genes)):
#     #can test fold diff and p-val for t-test
#     #present as volcano, dot size by present
#  #   print(genes[i],present[i],averages[i],totals[i])
#     if genes[i] in sigs:
#         plt.scatter(present[i],averages[i], c='red')
#         #plt.text(present[i],averages[i], genes[i])
#     else:
#         plt.scatter(present[i], averages[i], c='gray')
# plt.xlabel("Present")
# plt.ylabel("Averages")
# plt.title("Significant")
# #plt.scatter(t, p)
# plt.show()
# #%%
# """extra if you want it. Find how successful 
# is prediction based on any kind of selected data
# first, define which data to use and give indexes of ages
# """
# #run random forest
# X=np.array(selected).T[1:] #was data
# genes=np.array(selected).T[0] 
# y=age
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.datasets import make_classification
# clf = RandomForestClassifier(max_depth=2, random_state=0)
# clf.fit(X, y)
# print("success of prediction, score:", clf.score(X, y))
# imp=clf.feature_importances_
# sort=sorted(zip(imp, genes), reverse=True)
# print("Random forest age predictors")
# for i in range(len(sort))[:10]:
#     if sort[i][0]>0:
#         print(i, sort[i])
