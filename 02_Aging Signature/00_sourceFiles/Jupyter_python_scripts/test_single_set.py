"""
script test_single_set.py

derived from DGSE_sets_to_sign_scores
tests any other file compared to signature set.
any other file is added last to the collection of other sets
"""
import matplotlib.pyplot as plt, scipy.stats as ss, numpy as np

"input file with precalculated signature gene list"
f_in=open("Aging_List_REAN.txt","r")
lines=f_in.readlines()
f_in.close()
items=[line.rstrip().split("\t") for line in lines[1:]]
names, freqs =[x[0] for x in items], [int(x[1]) for x in items]
sort=sorted(zip(freqs,names), reverse=True)
mymax=sort[0][0]
print("mymax and all list",mymax, len(sort))

"Split sort list on gene names and gene scores"
GeneScore, Genes =[x[0] for x in sort], [x[1] for x in sort]
print("GeneScore", len(GeneScore), GeneScore) #a.k.a. rank score
print("Genes", len(Genes), Genes)

print("make control points for DGSE plot")
controls=[]
for i in range(0,max(GeneScore)):
    controls.append(len([x for x in GeneScore if x>i]))
control_points=controls[::-1] #reverse order
print(control_points) #[1, 4, 7, 16, 22, 57, 103, 185, 344, 875, 4024]
print(control_points[-1])
#calculate cumulative % in control points
for i in range(len(control_points)):
    print(i, control_points[i], control_points[i]*100/4024)


"single file for dGSE plot"
mdir="/Users/leonid_bystrykh/Documents/Articles&Presentations/2019/aging_sign/aging_HSCs_raw/meta_art/"
mfile="Noda.csv"
#mdir=""
#mfile="Bersenev_re.csv"

#RScores,EnScores, deltas=[],[],[] #all fishers and chisquares series, rank scores and deltas

f_in=open(mdir+mfile,"r")
lines=f_in.readlines()
f_in.close()
names=[]
for line in lines[2:]:
    items=line.rstrip().replace(" ","").split(",")
    item=items[0].split("///")
    names.append(item[0])

shortlist,scores=[],[] #gene indexes and scores
for name in names:
    if name in Genes:
        shortlist.append(Genes.index(name))
        scores.append(GeneScore[Genes.index(name)])
    else:
        scores.append(0) #unlisted genes go to the score 1
#print("check genes and scores found")
#for i in range(mymax+1):
#    print(i, scores.count(i))
"count rank score" #add all scores if >1, subtract all scores==1
RankScore=0
for s in scores:
    if s>1:
        RankScore+=s
    else:
        RankScore -=1
print("RankScore", RankScore)
#"make cumulative curve"
#observed,expected,delta,fisher,chisquare=[],[],[],[],[]
#for i in control_points:
#       # found = [x for x in shortlist if x < i]
#        found=[x for x in names if x in Genes[:i]]
#        print(i, "found", len(found))
#        observe = (len(found))
#        observed.append(observe)
#        expect = len(shortlist) * i / len(Genes)
#        expected.append(expect)
#        EnScore=((observe-expect)*100/len(names))
#        f_fi = ss.fisher_exact([[expect, len(shortlist) - expect], [observe, len(shortlist) - observe]])
#        print(i, len(shortlist), expect, observe, EnScore, f_fi[1], sep="\t")




