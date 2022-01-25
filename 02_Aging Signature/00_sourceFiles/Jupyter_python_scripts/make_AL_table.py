"""
take DE lists filtered by gene in signature
check them by gene list
try PCA
"""

import matplotlib.pyplot as plt, numpy as np
from mpl_toolkits.mplot3d import axes3d

#function for reading files
def getitems(fil):
    sfil=open(fil, "r"); lines=sfil.readlines(); sfil.close()
    items=[x.rstrip().split("\t") for x  in lines[1:]]
    return(items)

#load list and make signature
items=getitems("Aging_List_REAN.txt")
ASgenes=[x[0] for x in items if int(x[1])>1]
print("len items, len AL", len(items), ASgenes)

#load individual files and collect data
mdir="/Users/leonid_bystrykh/Documents/Articles&Presentations/2020/AgingSign/data/rean/AL_filtered/"
files=["Bersenev_AL.csv","Chambers_AL.csv","Flach_AL.csv","Norddahl_AL.csv","Wahlestedt_AL.csv",
       "Grover_AL.csv","Kirshner_AL.csv","Kowalczyk_AL.csv","Mann_AL.csv",
    "Lazare_AL.csv","Maryanovich_AL.csv","Sun_AL.csv"]
logFCs=[]
for file in files:
    print(file)
    myitems=getitems(mdir+file)
    mygenes=[x[1] for x in myitems]
    myfreqs=[float(x[2])for x in myitems]
    #collect logFCs, redundant gene names in arrays!
    freqs=[]
    for i in range(len(ASgenes)):
        fre=[]
        for j in range(len(mygenes)):
            if ASgenes[i]==mygenes[j]:
                fre.append(myfreqs[j])
        freqs.append(fre)
    #adjust to single value per gene
    new_freqs=[]
    for fre in freqs:
        if len(fre)>1:
            new_freqs.append(sum(fre)/len(fre))
        if len(fre)==1:
            new_freqs.append(fre[0])
        if len(fre)==0:
            new_freqs.append(0)
    logFCs.append(new_freqs)


def save_to(file_out):
    #file_out="AL_table.csv"
    f_out=open(file_out, "w")
    f_out.writelines("Files\Genes,"+ str(ASgenes)+"\n")
    for i in range(len(files)):
        f_out.writelines(files[i]+","+str(logFCs[i])[1:-1]+"\n")
    f_out.close()
def save_transposed():
    file_out="AL_transposed.csv"
    table=np.transpose(logFCs)
    print("table", table)
    f_out=open(file_out, "w")
    f_out.writelines("Files\Genes,"+ str(files)+"\n")
    for i in range(len(ASgenes)):
        f_out.writelines(ASgenes[i]+","+str(list(table[i]))[1:-1]+"\n")
    f_out.close()
save_transposed()

from sklearn.decomposition import PCA
fig = plt.figure()
ax = fig.add_subplot(111)#, projection='3d')

X = np.array(np.transpose(logFCs))
pca = PCA(n_components=3)
pca.fit(X)
pc1=pca.components_[0]
print("length PC", len(pc1))
pc2=pca.components_[1]
pc3=pca.components_[2]
colorz=["orange","orange","orange","orange","orange","brown","brown","brown","brown","red","red","red"]
#ax.scatter(pc1,pc2,pc3, c=colorz)
#plt.show()
#print(pca.components_)
print("explained", pca.explained_variance_ratio_)

for i in range(len(pc1)): #plot each point + it's index as text above
    ax.scatter(pc1[i],pc2[i],color=colorz[i])
    ax.text(pc1[i],pc2[i],  files[i][:10])
ax.set_xlabel('PC1 '+str(round(pca.explained_variance_ratio_[0],3)))
ax.set_ylabel('PC2 '+str(round(pca.explained_variance_ratio_[1],3)))
plt.show()

#for i in range(len(pc1)): #plot each point + it's index as text above
#    ax.scatter(pc1[i],pc2[i],pc3[i],color=colorz[i])
#    ax.text(pc1[i],pc2[i],pc3[i],  files[i][:10])
#ax.set_xlabel('PC1')
#ax.set_ylabel('PC2')
#ax.set_zlabel('PC3')
#plt.show()

