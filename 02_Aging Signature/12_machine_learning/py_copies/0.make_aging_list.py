#!/usr/bin/env python
# coding: utf-8

# Probably the first script to use.
# You will need series of files with lists of differentially expressed genes.
# This script will compile a unified list of all genes from all lists and give a citation frequency of each gene. By our definition (used through the paper) all genes with frequency above 3 comprise the Aging Signature list. 
# 
# Make nesessary preparations for the analysis of the aging data. You will need:
# 1. Recent python 3+ version
# 2. Installed jupyter notebook (via anaconda or standalone)
# 3. Source files. The files used in the paper can be found next to this script. They are all uniformly organised as tab-separated text files. Each of them contains two lines of headings. Gene names are in the second column, fold changes are in the third column. All other data are not used in current scripts.
# Note that python counts from 0, therefore the second column is *col[1]*, and the 3rd column is *col[2]* 

# Import necessary packages

# In[1]:


import matplotlib.pyplot as plt, numpy as np, math, scipy.stats as stats 
print("imported")


# Define directory where you keep your source files: "mdir"
# List all files you will import for analysis: "sets"
# Make shortlists from all signature files: "setnames"
# The block below is for META-files, meaning they are derived from the supplementary files provided by authors

# In[2]:


#meta data from meta directory
sets=["Beerman_2013.txt","Chambers_2007.txt","Flach_2014.txt","Grover_2016.txt",
    "Kirshner_2017.txt","Kowalczyk_2015.txt","Noda_2009.txt","Quere_2014.txt",
    "Rossi_2005.txt","Sun_2014.txt","Wahlestedt_2013.txt"]
setnames=[x[:-4] for x in sets]
"add  other list to the list of other files"
"extra file for DGSE plot and comparison to others"

mdir="/Users/leonid/Documents/Publications_Presentations/2020/AgingSign/data/meta/"


# The block below is for reanalysed files.  The structure and variable names are the same as above, therefore you can use either of those block, but not both at the same time (this is not recommended, because some of the files will be partial duplicates of itself and cause statistical bias).

# In[2]:


#reanalysed files from rean directory
sets=["Bersenev_GSE39553.csv","Chambers_GSE6503.csv",
    "Flach_GSE48893.csv","Grover_GSE70657.csv",
    "Kirshner_GSE87631.csv","Kowalczyk_GSE59114.csv",
    "Lazare_GSE128050.csv","Mann_GSE1004426.csv",
      "Maryanovich_GSE109546.csv","Norddahl_GSE27686.csv",
      "Sun_GSE47817.csv","Wahlestedt_GSE44923.csv"]
setnames=[x[:-4] for x in sets]
"add  other list to the list of other files"
"extra file for DGSE plot and comparison to others"

mdir="/Users/leonid/Documents/Publications_Presentations/2020/AgingSign/data/rean/"


# open every file and collect all names of genes
# names is set of sets, 
# fcs are fold changes per gene
# allg (all genes) and 
# allf (all found times) are united set with frequencies

# In[3]:


names,fcs,allg,allf=[],[],[],[]
print("file\t\t\tlength\tpositive\tnegative")
for file in sets:
    f_in=open(mdir+ file,'r')
    lines=f_in.readlines()
    na,fc=[],[] #list of genes (na) and frequencies (fc) per set
    for line in lines[2:]: #skip first two lines
        items= line.rstrip().replace("'","").replace(" ","").split('\t')
        it=items[1].split('///') #for multiple gene names take the 1st
        if it[0] not in na:
            na.append(it[0]) #0 stands for the 1st column, otherwise change
            fc.append(float(items[2])) #1 stands for the second column
    #count all genes and their frequencies per publication
    for n in na:
        if n not in allg:
            allg.append(n)
            allf.append(1)
        else:
            allf[allg.index(n)]+=1
    names.append(na) #list of lists of gene names
    fcs.append(fc)  #list of lists of fold changes
    pos=[x for x in fc if x>0]
    neg=[x for x in fc if x<0]
    print(file, len(fc), len(pos), len(neg), sep="\t")
    f_in.close()


# All requred data are collected.
# We can check every set name, lenth of every gene set, print first 10 names from each set
# We also check for possible damage in gene names.
# Some tables, when open/saved in excel loose original gene names, which became a dates, like "Mar-2", or "Sep-5". so we collect those if they are (or any other characters)

# In[4]:


print("summary")
for i in range(len(sets)):
        print(sets[i],len(names[i]),names[i][:5],"...")
        mar=[x for x in names[i] if x[:4]=="Sep-"]
        if len(mar)>0:
            print("found Mar", sets[i], mar)
    


# Next, we collect all frequencies for all genes. 

# In[5]:


allfcs=[] #all fold changes
for gene in allg:
    folds=[]
    for i in range(len(names)):
        for j in range(len(names[i])):
            if gene==names[i][j]:
                folds.append(fcs[i][j])
    allfcs.append(folds)
print("ready to make the AgingList file")


# Now we can save the data in the Aging List file.
# Choose one of two blocks below, depending on what set of data you have used. We will save average of all fold changes.
# We can also save all fold changes as a list in a separate column, if we want. To do so just add allfcs[i] as separate column in f_out.writelines() line
# 
# block below is for META data, and the list will be saved into "rean_gene_fcs_META.txt".

# In[6]:


f_out=open("Aging_List_META.txt","w")
f_out.writelines("Gene\tFreq_group\tLog2FC\n")

for i in range(len(allg)):
   # print(allg[i],allf[i], allfcs[i])
    f_out.writelines(allg[i]+"\t"+str(allf[i])+"\t"+str(sum(allfcs[i])/len(allfcs[i]))+"\n")
f_out.close()


# block below is for reanalysed data, and the list will be saved into "rean_gene_fcs_REAN.txt".

# In[9]:


f_out=open("Aging_List_REAN.txt","w")
f_out.writelines("Gene\tFreq_group\tLog2FC\n")

for i in range(len(allg)):
   # print(allg[i],allf[i], allfcs[i])
    f_out.writelines(allg[i]+"\t"+str(allf[i])+"\t"+str(sum(allfcs[i])/len(allfcs[i]))+"\n")
f_out.close()
print("AgingList saved")


# In principle, the Aging List is the most important part of the scipt. You can stop here and proceed with the next script for further study and use of the aging signature. 
# Below are some additional files which you can use if you want it.
# Here is a cross-comparison of all files used for the list.

# In[7]:


print("table of intersects\t set_size\t Intersects_tab\t all_intersected\t nr_intersected\t fraction_intersected")
s,v,u,dists=[],[],[],[] # u for unique genes
for i in range(len(names)):
    val,uni=[],[]
    "define pos, neg"
    for j in range(len(names)):
        ovl=[x for x in names[i] if x in names[j]]
        if i!=j:
            val.append(len(ovl))
            for o in ovl:
                if o not in uni:
                    uni.append(o)
        else:
            self=len(ovl) #dirty trick for own size
            val.append(0)
    u.append(uni)
    s.append(self)
    dists.append(val)
    v.append(len(uni)) #length unique list
    print(sets[i],self,val,sum(val),len(uni), len(uni)/self, sep="\t")
    


# The table above can be presented as a heatmap below

# In[8]:


#heatmap figure for similarities between gene sets
plt.imshow(dists)
plt.yticks(range(len(setnames)),setnames)
plt.colorbar()
#plt.show()


# Scatter plot can be skipped
# At the end some summary table for each file used for the Aging List

# In[9]:


#plot of gene frequencies
#plt.plot(freqs,counts, 'o-')
#plt.scatter(freqs,counts)
#plt.show()

print("make groups of genes by their frequencies of occurrence")
print("per set")
print("set\ttotal\t>3\t relative") #\t p_hypergeometric
rel_dists=[]
for val in range(len(sets)):
    total,sig=len(names[val]),0
    for name in names[val]:
        if name in allg:
            score=allf[allg.index(name)]
            if score>3:
                sig+=1
    rel_dists.append(float(sig)/float(total))
    print(sets[val], total, sig, float(sig)/float(total), sep="\t") #myp,


# We can draw those results in different ways. Here is the script to make a "dart-plot"
# of data sets similarity degree to the aging signature. To do so we need to use Set names, their sizes, and relative distances to the Aging Signature. Here we use the last column from the table above, namely relative values (number of genes found in signature divided by total number of genes.

# In[13]:


sizes=[len(x) for x in names]
for i in range(len(sizes)):
    print(sets[i], sizes[i],rel_dists[i], sep="\t")


# if sets, sizes and rel_dists are OK, we can do a little geometric exersise and draw the "dart-Plot". First, we make some scale adjustments to the data. It is needed to adjust it to the figure size

# In[14]:


u,v=[],[]
size=[float(x)/10 for x in sizes]
dist=[(1-float(x))*10 for x in rel_dists]


# next we use our geometric knowledge and find positions for the nodes we will draw

# In[15]:


setlength=len(dist)
for i in range(setlength):
    u.append(round(math.sin(i*2*np.pi/setlength)*float(dist[i]),4))
    v.append(round(math.cos(i*2*np.pi/setlength)*float(dist[i]),4))
print(u)
print(v)


# Now we can draw our nodes and edges. Note that colors and line thickness can be changed too

# In[16]:


from matplotlib.patches import Ellipse
import numpy.random as rnd

#my struggles to control the figure size
#plt.figure(figsize=(15,5))
fig, axes= plt.subplots(nrows=1, ncols=1,figsize=(15,10))
#fig = plt.figure(0)
ax = fig.add_subplot(111)#, aspect='equal')
"draw ellipses"
for i in range(len(u)):
    ells = [Ellipse(xy=(u[i],v[i]), width=size[i]/30, height=size[i]/30,
                    facecolor=rnd.rand(3),alpha=0.9, edgecolor="black")]
    for e in ells:
        ax.add_artist(e)
"add lines"
for i in range(len(v)):
    plt.plot([0,u[i]],[0,v[i]],"--", color="black")
    plt.text(u[i],v[i],sets[i][:-7])
"add central hub"
ell=Ellipse(xy=(0,0), width=1, height=1,facecolor="red",alpha=1, edgecolor="black")
ax.add_artist(ell)

ax.set_xlim(min(u+v)-20, max(u+v)+20)
ax.set_ylim(min(u+v)-10, max(u+v)+10)

#ax.set_aspect(0.9) #shape of the graph rectangle
plt.gca().axes.get_yaxis().set_visible(False)
plt.gca().axes.get_xaxis().set_visible(False)


# Depending on your browser you might see the same figures differently.  Details can be found in stackoverflow: https://stackoverflow.com/questions/19410042/how-to-make-ipython-notebook-matplotlib-plot-inline
# Further adjustments are probably better to make in regular IDEs, like Spyder3 or PyCharm CE. 
# 
