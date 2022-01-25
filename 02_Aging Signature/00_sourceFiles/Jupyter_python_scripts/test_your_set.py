"""
derived from DGSE_sets_to_sign_scores
tests any other file compared to signature set.
any other file is added last to the collection of other sets
"""
import matplotlib.pyplot as plt, scipy.stats as ss, numpy as np

"input file with precalculated signature gene list"
f_in=open("rean_gene_fcs_art.txt","r")
lines=f_in.readlines()
f_in.close()
names, freqs=[],[]
for line in lines[1:]:
    items=line.rstrip().split("\t")
    names.append(items[0])
    freqs.append(int(items[1]))
print("names", len(names))
print("freqs", len(freqs))
sort=sorted(zip(freqs,names), reverse=True)
print(sort)
mymax=sort[0][0]
print("mymax",mymax)

"Split sort list on gene names and gene scores"
GeneScore,Genes=[],[]
for item in sort:
    GeneScore.append(item[0])
    Genes.append(item[1])
#for i in range(mymax+1):
 #   print(i, GeneScore.count(i))

print("GeneScore", len(GeneScore), GeneScore) #a.k.a. rank score
print("Genes", len(Genes), Genes)

print("make control points for DGSE plot")
control_points=[]
for i in range(0,len(GeneScore)-1):
    if GeneScore[i]!=GeneScore[i+1]:
        control_points.append(i+1)
        print(sort[i])
control_points.append(len(GeneScore)) #add last control point
print("Control points", control_points)
print("Control point 7", control_points[7])

"make shortlists from all signature files"
files=['Bersenev_re.csv','Chambers_re_em.csv','Flach_rmdup_re.csv','Grover_re.csv',
    'Kirshner_re.csv','Kowalczyk_re.csv','Mann_re.csv','Maryanovich_re.csv',
    'Norddahl_rmdup_re.csv','Seka_re.csv','Sun_re.csv','Wahlestedt_rmdup_re.csv'
      ]
setnames=['Bersenev','Chambers','Flach','Grover',
    'Kirshner','Kowalczyk','Mann',"Maryanovich",'Norddahl','Lazare',
    'Sun','Wahlestedt'
      ]
"add  other list to the list of other files"
"extra file for DGSE plot and comparison to others"
mdir="/Users/leonid_bystrykh/Documents/Articles&Presentations/2019/aging_sign/aging_HSCs_raw/meta_art/"
mfile="Noda.csv"
#mdir=""
#mfile="signature_list_max.csv"
files.append(mdir+mfile)
setnames.append(mfile[:-4])

"run through all files and estimate scores"
print("file", "observed","expected","delta","enrichment%", "fisher","chisquare", sep="\t")
fis,chis,RScores,EnScores, deltas, set_len=[],[],[],[],[] ,[]#all fishers and chisquares series, rank scores and deltas
for file in files:
    f_in=open(file,"r")
    lines=f_in.readlines()
    f_in.close()
    names=[]
    for line in lines[2:]:
        items=line.rstrip().replace(" ","").split(",")
        item=items[0].split("///")
        names.append(item[0])
    shortlist,scores=[],[] #gene indexes and scores
    set_len.append(len(names))
    #calculate rank score
    for name in names:
        if name in Genes:
            shortlist.append(Genes.index(name))
            scores.append(GeneScore[Genes.index(name)])
          #  if file=="signature_list_max.csv":
           #     print("max score list", name, GeneScore[Genes.index(name)])
        else:
            scores.append(1) #unlisted genes go to the score 1
         #   if file=="signature_list_max.csv":
          #      print(name, " not in the list")
    #print("list lengthes", file, len(names), len(shortlist))
    "count rank score" #add all scores if >1, subtract all scores==1
    RankScore=0
    for s in scores:
        if s>1:
            RankScore+=s
        else:
            RankScore -=s
    RScores.append(RankScore) #needed for pin-plot
    "make cumulative curve"
    observed,expected,delta,fisher,chisquare=[],[],[],[],[]
    for i in control_points:
            found = [x for x in shortlist if x < i]
            observe = (len(found))
            observed.append(observe)
            expect = len(shortlist) * i / len(Genes)
            expected.append(expect)
            delta.append(observe - expect)
            f = ss.chisquare([observe, len(shortlist) - observe], f_exp=[expect, len(shortlist) - expect])
            chisquare.append(f[1])
            f = ss.fisher_exact([[expect, len(shortlist) - expect], [observe, len(shortlist) - observe]])
            fisher.append(f[1])
    EnScores.append(delta[7]*100/len(names)) #Enrichment scores at 5%
    print(file, observed[7],expected[7],delta[7],EnScores[-1],fisher[7],chisquare[7])
    fis.append(fisher)
    chis.append(chisquare)
    deltas.append(delta)
print("Sets, gene_size,  RankScore, EnScore")
for i in range(len(files)):
    print(setnames[i],set_len[i], RScores[i], EnScores[i])

plt.figure(figsize=(16,7))
"make DGSE plot for the last set"
plt.subplot(121)
plt.title(setnames[-1])
delt=deltas[-1] #absolute bias in gene numbers, take the last series
plt.plot(control_points, delt)
plt.ylabel("Normalized enrichment, %")
#red ticks
for i in range(len(shortlist)):
        x = [shortlist[i], shortlist[i]]
        y = [0, -5]
        plt.plot(x, y, c="red")

"make 1D pinplot for Rank scores"
plt.subplot(222)
data, scores,highs=fis,[],[]
scores=RScores
for i in range(len(scores)):
    if i==len(scores)-1: #if the last one
        height = 10
        plt.plot([scores[i],scores[i]], [0,height],c="black", alpha=0.5) #line
        plt.text(scores[i],height,setnames[i]) #text
        plt.scatter(scores[i],height, s=100, c="red",alpha=1) #dot
    else:
        height=5
        plt.plot([scores[i], scores[i]], [0, height], c="red", alpha=0.5)  # line
        plt.text(scores[i], i/2, setnames[i])  # text
        plt.scatter(scores[i], height, s=100, c="black", alpha=1)  # dot
plt.xlabel("Rank scores")
plt.yticks([],[])

"pinplot for Enrichment scores"
plt.subplot(224)
data, scores,highs=fis,[],[]
scores=EnScores
for i in range(len(scores)):
    if i==len(scores)-1: #if the last one
        height = 10
        plt.plot([scores[i],scores[i]], [0,height],c="black", alpha=0.5) #line
        plt.text(scores[i],height,setnames[i]) #text
        plt.scatter(scores[i],height, s=100, c="red",alpha=1) #dot
    else:
        height=5
        plt.plot([scores[i], scores[i]], [0, height], c="red", alpha=0.5)  # line
        plt.text(scores[i], i/2, setnames[i])  # text
        plt.scatter(scores[i], height, s=100, c="black", alpha=1)  # dot
plt.xlabel("Enrichment scores")

#plt.savefig(file[:-4]+".pdf")
plt.yticks([],[])
plt.show()


