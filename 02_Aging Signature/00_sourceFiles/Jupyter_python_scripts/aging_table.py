"""
compare lists overlaps and distances

"""
import matplotlib.pyplot as plt, numpy as np, math, scipy.stats as stats #, skbio.diversity.alpha as sda
"load signature"
md="/Users/leonid_bystrykh/Documents/Articles&Presentations/2020/AgingSign/data/aging_py/"
asfl="Aging_List_REAN.txt"
f_in=open(md+asfl, "r")
lines=f_in.readlines()
genes, freqs=[],[]
for line in lines[1:]:
    items=line.rstrip().split("\t")
    genes.append(items[0])
    freqs.append(items[1])
print("stats of frequencies")
for item in set(freqs):
    print(item, freqs.count(item))
#META data
M_sets=["Beerman_2013.txt","Chambers_2007.txt","Flach_2014.txt","Grover_2016.txt",
    "Kirshner_2017.txt","Kowalczyk_2015.txt","Noda_2009.txt","Quere_2014.txt",
    "Rossi_2005.txt","Sun_2014.txt","Wahlestedt_2013.txt"]
M_dir="/Users/leonid_bystrykh/Documents/Articles&Presentations/2020/AgingSign/data/meta/"

R_dir="/Users/leonid_bystrykh/Documents/Articles&Presentations/2020/AgingSign/data/rean/"
"make shortlists from all signature files"
R_sets=["Bersenev_GSE39553n_DE_YOn.csv",
    "Chambers_GSE6503_YOn.csv",
    "Flach_GSE48893_YOn.csv",
    "Grover_GSE70657_DE_OY.csv",
    "Kirshner_GSE87631_DE_OY.csv",
    "Kowalczyk_GSE59114_DE_OY.csv",
    "Lazare_DE_YO.csv",
    "Mann_GSE1004426_DE_OY.csv",
    "Maryanovich_GSE109546_DE_YO.csv",
    "Norddahl_GSE27686_DE_YOn.csv",
    "Sun_GSE47817_DE_YO.csv",
    "Wahlestedt_GSE44923_DE_YOn.csv"
      ]
sets, mdir=R_sets, R_dir

setns=[x.split("_") for x in sets]
setnames=[x[0] for x in setns]
"add  other list to the list of other files"

"names is set of sets, allg (all genes) and allf (all found times) are united set with frequencies"
names,fcs,allg,allf=[],[],[],[]
for set in sets:
    f_in=open(mdir+ set,'r')
    lines=f_in.readlines()
    na,fc=[],[] #list of genes and frequencies per set
    for line in lines[2:]: #skip first two lines
        if "Gm4604" in line:
            print("\nfound name",set, line,"\n")
        items= line.rstrip().replace("'","").replace(" ","").split('\t')
        it=items[0].split('///') #for multiple gene names take the 1st
        if it[0] not in na:
            na.append(it[0])
            fc.append(float(items[2]))

    #count all genes and their frequencies per publication
    for n in na:
        if n not in allg:
            allg.append(n)
            allf.append(1)
        else:
            allf[allg.index(n)]+=1
    names.append(na) #list of lists of gene names
    fcs.append(fc)  #list of lists of fold changes
    f_in.close()

#f_out=open("rean_gene_fcs_META.txt","w")
#f_out.writelines("Gene\tFreq_group\tLog2FC\n")
#allfcs=[] #all fold changes
#for gene in allg:
#    folds=[]
#    for i in range(len(names)):
#        for j in range(len(names[i])):
#            if gene==names[i][j]:
#                folds.append(fcs[i][j])
#    allfcs.append(folds)
#for i in range(len(allg)):
#   # print(allg[i],allf[i], allfcs[i])
#    f_out.writelines(allg[i]+"\t"+str(allf[i])+"\t"+str(sum(allfcs[i])/len(allfcs[i]))+"\n")
#f_out.close()

print("summary")
for i in range(len(sets)):
    print(sets[i],len(names[i]),names[i])
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
#heatmap figure for similarities between gene sets
#plt.imshow(dists)
#plt.yticks(range(len(setnames)),setnames)
#plt.colorbar()
#plt.show()
print("make groups of genes by their frequencies of occurrence")

#plot of gene frequencies
#plt.plot(freqs,counts, 'o-')
#plt.scatter(freqs,counts)
#plt.show()

print("per set")
print("set\ttotal\t>3\t relative") #\t p_hypergeometric
for val in range(len(sets)):
    total,sig=len(names[val]),0
    for name in names[val]:
        if name in allg:
            score=allf[allg.index(name)]
            if score>3:
                sig+=1
    print(sets[val], total, sig, float(sig)/float(total), sep="\t") #myp,
