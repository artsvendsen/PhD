"""
derived from test_your_set.py
tests only perfect top5 signature list and random list from all reported genes

"""
import matplotlib.pyplot as plt, random

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

"find top 5% range and collect top 5% gene names"
print("Freq groups")
listlen, top5 =0,[x for x in freqs if x>2]
for i in range(len(freqs)):
    myselection=len([x for x in freqs if x>i])
    print(i+1, myselection, 100*myselection/len(freqs))

"make random subset i times"
allscores,allenscores=[],[]
for i in range(100):
    mylen=random.choice(range(100,500))
    my_selection=[]
    while len(my_selection)<mylen:
        my_selection.append(random.choice(freqs))
    rank_score=0
    for name in my_selection:
        if name>1:
            rank_score+=name
        else:
            rank_score -= name
    "calculate enrichment score"
  #  print("i, observed, expected, EnScore, RankScore")
    expect=len(my_selection)*len(top5)/len(freqs)
    observed=len([x for x in my_selection if x >2])
    allscores.append(rank_score)
    en_score=(observed-expect)*100/len(my_selection)
    allenscores.append(en_score)
    print(i, mylen, expect, observed, en_score,rank_score, sep="\t")
plt.subplot(211)
plt.hist(allscores)
plt.xlabel("Rank_scores")
plt.subplot(212)
plt.hist(allenscores)
plt.xlabel("Enrich_scores")
plt.show()
