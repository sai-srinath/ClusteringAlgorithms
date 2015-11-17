import csv
import numpy as np
from sklearn.decomposition import PCA 
import pylab as pl
import sys
path=sys.argv[1]
path2=sys.argv[2]
reader = csv.reader(open(path),delimiter="\t")
l = list(reader)
for i in l:
    del i[0]
    del i[1]
    
reader2= reader = csv.reader(open(path2),delimiter="\t")
target=list(reader2)   

l=np.array(l)
l=l.astype(np.float)
target=np.array(target)
target=target.astype(np.float)
print l.shape
print target.shape

pca=PCA(n_components=2)
pca.fit(l)
l_reduced=pca.transform(l)

pl.scatter(l_reduced[:,0],l_reduced[:,1],s=100,c=target)  
pl.show() 

