#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 16:18:41 2021
Sankey did not work,
I made this kind of a transition plot instead
@author: leonid
"""

import pandas as pd, os, matplotlib.pyplot as plt
#%%
os.chdir("/Users/Art/Drive/PhD/Experiments/Aging Signature/000_FINAL/12_machine_learning/selected_AS/")
os.getcwd()
df1=pd.read_csv("Best_ADAboost.tsv", index_col=0, sep="\t")
df2=pd.read_csv("Best_RNDforest.tsv", index_col=0, sep="\t")
print("index file 1", df1.index,"\nIndex file2", df2.index)
lim=50
words1=df1.index[:lim]
words2=df2.index[:lim]
#%% generate y coords for lines and arrows
Y=[]
set1=[-1]*max(len(words1),len(words2))
set2=[-1]*max(len(words1),len(words2))
for i in range(len(words1)):
    for j in range(len(words2)):
        if words1[i]==words2[j]:
            Y.append([i,j])

#%%
#generate start and stop positions
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['svg.fonttype'] = 'none'
for i in range(len(Y)):
    if min(Y[i])>=0:
        y=Y[i]
        x=[1,2]
        colr="gray"
        if y[0]>y[1]:
            colr="red"
        if y[0]<y[1]:
            colr="blue"
        plt.plot(x,y, c=colr, linewidth=2)
     #   plt.arrow(x[0],y[0],x[1]-x[0],y[1]-y[0],
      #            head_width=0.1, color=colr, linewidth=3)
     #   plt.text(x[0]-0.4,y[0],words[i], fontsize=15, color=colr)
      #  plt.text(x[1]+0.3,y[1], words[i], fontsize=15, color=colr)
plt.gca().invert_yaxis()
# words will be positioned according to their indexes
for i in range(len(words1)):
    plt.text(0.8,i,words1[i])
for i in range(len(words2)):
    plt.text(2,i,words2[i])
#plt.arrow(1,1,1,1, head_width=1)
plt.xticks([],[])
plt.yticks([],[])
plt.xlim(0,3)
plt.ylim(lim,-1)
plt.title("Best ADAboost vs RNDforest")
plt.show()
#plt.savefig("test.svg", format="svg")
# %%
