import os, sys, numpy, math
from MolKit2 import Read
from mglutil.math.kmeansClustering import kmeans, Point

mol = Read(sys.argv[1])
#numClusters = int(sys.argv[2])
basename = os.path.splitext(sys.argv[1])[0]

coords = mol._ag._coords
# nc is the number of models, na number of atoms
nc, na, dum = coords.shape

cg = [numpy.sum(coords[i], 0)/na for i in range(nc)]

from DBSCAN import dbscan
eps = 2.0
minpts = 2
pointlabel,cl = dbscan(cg, eps, minpts)
nbCluster = max(pointlabel)+1
clusters = []

for n in range(nbCluster):
    clusters.append([])

for posecgind, cln in enumerate(pointlabel):
    clusters[cln].append(posecgind+1)

for cl in clusters:
    print(len(cl), numpy.array(cl))

import prody
from mglutil.util.io import StringIO
for n, cl in enumerate(clusters):
    if not os.path.exists('%s_clusters'%basename):
        os.mkdir('%s_clusters'%basename)
    f = open('%s_cl%03d.pdb'%(basename, n+1), 'w')
    for c in cl:
        if not os.path.exists('%s_clusters/%d'%(basename, n)):
            os.mkdir('%s_clusters/%d'%(basename,n))
        prody.writePDB('%s_clusters/%d/%s_pose%04d'%(basename,n,basename,c),
                       mol._ag, csets=[c-1])
        buf = StringIO()
        prody.writePDBStream(buf, mol._ag, csets=[c-1])
        f.write('MODEL   %6d\n'%c)
        lines = buf.readlines()[1:]
        [f.write('%s'%l) for l in lines]
        f.write('ENDMDL\n')
    f.close()
    
## f = open('wr6cg.npy', 'w')
## numpy.save(f, cg)
## f.close()

## pts = [Point(x) for x in cg]

## clusters = kmeans(pts, numClusters, 2.0, initial=None)


## clusterInds = []
## for cl in clusters:
##     inds = []
##     for pt in cl.points:
##         x = pt.coords
##         inds.append(lookup['%.3f %.3f %.3f'%(x[0],x[1],x[2])])
##     clusterInds.append(inds)

## for ci in clusterInds: print("len: %d %s"%(len(ci), ci))
