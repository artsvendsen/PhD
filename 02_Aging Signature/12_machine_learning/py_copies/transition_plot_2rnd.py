#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 17:06:22 2021
kind of replacement of sankey
https://blog.ouseful.info/2017/11/28/quick-round-up-visualising-flows-using-network-and-sankey-diagrams-in-python-and-r/
better known as transition matrix
@author: leonid
"""
#%%
import matplotlib.pyplot as plt, random

#take one list, reshuffle and make the second one
words=["one","two", "three","four","пять","6"]
words1,words2=[],[]
#while len(words2)<5:
for i in range(10):
    words1.append(random.choice(words))
    words2.append(random.choice(words))
print("1st set of words", words1)
print("2nd set of words", words2)
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
for i in range(len(Y)):
    if min(Y[i])>=0:
        y=Y[i]
        x=[1,2]
        colr="gray"
        if y[0]>y[1]:
            colr="red"
        if y[0]<y[1]:
            colr="blue"
        plt.plot(x,y, c=colr, linewidth=3)
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
plt.ylim(len(max(words1,words2)),-1)
plt.show()
# %%
