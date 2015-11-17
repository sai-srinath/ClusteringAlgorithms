import csv
import numpy as np
from math import sqrt
from scipy.spatial import distance
from scipy.cluster.hierarchy import dendrogram, linkage
import math
import matplotlib.pyplot as plt
from sklearn.metrics import jaccard_similarity_score
import pylab as pl
from sklearn.decomposition import PCA 
import sys

    
path = sys.argv[1]
numberofclusters=sys.argv[2]
reader = csv.reader(open(path),delimiter="\t")
l = list(reader)
#for pca
#pcal=l
#for dendogram
X=l

#for Jaccard
noc=10
true_labels=[]
for i in range(noc):    
    true_labels.append([])
    
pointmap = {}
for v in l:
    #removing the true label
    if v[1]!='-1':
        true_labels[int(v[1])-1].append(v[0])
    
    del v[1]
    #storing the point ID
    val=v[0]
    #removing point ID
    del v[0]
    pointmap[val]=np.array(v,dtype=float)


#distance.euclidean(pointmap['478'],pointmap['479'])   


for v in pointmap:
    pointmap[v]=np.asmatrix(pointmap[v])


#make distance matrix
size=len(pointmap)
distmatrix = np.zeros((size,size),dtype=float)

indexmap={}
index=0;
for k in pointmap:
    indexmap[index]=k
    index=index+1
#initializing distance matrix
for i in range(size):
    for j in range(size):
        ID1=indexmap[i]
        ID2=indexmap[j]
        distmatrix[i][j]=distance.euclidean(pointmap[ID1],pointmap[ID2])
        #print i,j,distmatrix[i][j]

#must be iterated
emptyct=0
while emptyct!=(size-int(numberofclusters)-1):
    emptyct=0
    for k in indexmap:
        if indexmap[k]=='empty':
            emptyct=emptyct+1
    #print "number of emptyct is:",emptyct
    
                    
    minimum=distmatrix[distmatrix>0].min()
    k=np.where(distmatrix==minimum)[0]
    
    
    for i in range(size):
        #row-wise comparison
        distmatrix[k[0]][i]=min(distmatrix[k[0]][i],distmatrix[k[1]][i])
        distmatrix[k[1]][i]=0
        #column-wise comparison
        distmatrix[i][k[0]]=min(distmatrix[i][k[0]],distmatrix[i][k[1]])
        distmatrix[i][k[1]]=0
    
    #merging the IDs 
    #indexmap[k[0]]=indexmap[k[0]]+","+indexmap[k[1]]
    #indexmap[k[1]]="empty"
    
    
    #indexmap[k[0]].append(val)
    print "merging",indexmap[k[0]],"and ",indexmap[k[1]]
    indexmap[k[0]]=indexmap[k[0]]+","+indexmap[k[1]]
    
    
    indexmap[k[1]]="empty"
    #print "after merging",indexmap[k[0]]
    
    #print indexmap[k[0]].count(',')+1
    
#dendogram
for i in X:
    del X[0]
    del X[1]

   
Z=linkage(X,'single')
plt.figure()
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Index')
plt.ylabel('Euclidean DIstance')
dendrogram(Z,leaf_rotation=90.,leaf_font_size=8.,)
plt.show()

#finding Jaccard coefficient
cluster_results=[]
for i in range(noc):    
    cluster_results.append([])
#print "printing indexmap"
ct=0
for i in range(size):
    if indexmap[i]!='empty':
        cluster_results[ct].append(indexmap[i])
        ct+=1
#print len(cluster_results)

# 
newlist=[]
for i in range(len(cluster_results)):
    newlist.append([])
for i in range(len(cluster_results)):    
    a=cluster_results[i]
    tmp=[]
    for j in a:
        for v in j.split(','):
            tmp.append(v)
    newlist[i]=tmp   

#newlist stores the results of clusters                           
    
true = np.zeros((size,size),dtype=int)
result = np.zeros((size,size),dtype=int) 
#initializing true label matrix
for i in range(size):
    for labels in true_labels:
        if str(i+1) in labels:
            for l in labels:
                true[i][int(l)-1]=1  

#initializing final cluster matrix
for i in range(size):
    for labels in newlist:
        if str(i+1) in labels:
            for l in labels:
                result[i][int(l)-1]=1                    
   
#print result   

score=jaccard_similarity_score(true,result)  
print "Jaccard coefficient is"
print score  

#finding PCA


#target values:
target=np.zeros((size),dtype=np.float)
#for i in range(size):
#    target.append([])
for i in range(1,size):
    for j in range(len(newlist)):
        if str(i) in newlist[j]:
            #print "yes",str(i),"is there"
            target[i]=j

reader2 = csv.reader(open(path),delimiter="\t")
l2 = list(reader2)

for i in l2:
    del i[0]
    del i[1]
    
pcal=l2

pl.figure()
pcal=np.array(pcal)
pcal=pcal.astype(np.float)  
#print pcal.shape
#print target.shape  
pca=PCA(n_components=2)
pca.fit(pcal)
pcal_reduced=pca.transform(pcal)
pl.scatter(pcal_reduced[:,0],pcal_reduced[:,1],s=50,c=target)  
pl.show() 
#writing target to file
f=open('target_hierarchical_iyer_10.txt','w')
for i in target:
    s=int(i)+1
    f.write(str(s)+'\n')
f.close() 

  

    