#!/u/home/s/sanjaydm/.local/bin/python3.5
from math import *
from numpy import dot, array
import sys
fileCyl = "./cylinder_4260.vtk"
fileBump="./bump.vtk"
fileBumpGhost="./bumpGhost.vtk"
filenodesdat = "./gauss_nodes.dat"
fileconndat ="./gauss_conn.dat"
prefix=""
eps = float(sys.argv[1])
print("eps = %f" %eps)

def loadNodesFromVTK(filename):
    fid = open(filename,'r')
    X=[];Y=[];Z=[]
    # Read header
    for i in range(4):
        fid.readline()
    # The fifth-line contains info about number of nodes
    [PTS,numNodes,TYPE] = (fid.readline()).split()
    numNodes = int(numNodes)
    # Read the whole file as a list of strings
    XYZstr = (fid.read()).split()
    cntr = 0
    while cntr < numNodes:
        X.append(float(XYZstr[3*cntr + 0]))
        Y.append(float(XYZstr[3*cntr + 1]))
        Z.append(float(XYZstr[3*cntr + 2]))
        cntr = cntr + 1
    fid.close()
    return [X,Y,Z]

def loadConnFromVTK(filename):
    fid = open(filename,'r')
    conn = [];
    # Read header
    for i in range(4):
        fid.readline()
    
    # The fifth-line contains info about number of nodes
    [PTS,numNodes,TYPE] = (fid.readline()).split()
    numNodes = int(numNodes)
    # Read the whole file as a list of strings
    dataStr = (fid.read()).split()
    # Read size of connectivity table info stored at index '3numNodes+1'
    numEle = dataStr[3*numNodes+1]
    numEle = int(numEle)
    # Note: 3N = POLYGON
    # Note: 3N+1 = int(numEle)
    # Note: 3N+2 = int(numEle*4)
    strt = 3*numNodes+3
    cntr = 0
    while cntr < numEle:
        three = int(dataStr[strt+4*cntr + 0])
        A = int(dataStr[strt+4*cntr + 1])
        B = int(dataStr[strt+4*cntr + 2])
        C = int(dataStr[strt+4*cntr + 3])
        conn.append([A, B, C])
        cntr = cntr + 1
    
    fid.close()
    return conn


def writeXYZToVTK(fout_vtk, X,Y,Z):
    if len(X) != len(Y) or len(X) != len(Z):
        print("Length mismatch")
        return -1
    numNodes = len(X)
    fout_vtk.write("# vtk DataFile Version 3.1\n")
    fout_vtk.write("vtk output\n")
    fout_vtk.write("ASCII\n")
    fout_vtk.write("DATASET POLYDATA\n")
    fout_vtk.write("POINTS "+str(numNodes)+" FLOAT\n")
    for i in range(numNodes):
        fout_vtk.write(str(X[i])+" "+str(Y[i])+" "+str(Z[i])+"\n")

def writeXYZToNodesDat(fout_vtk, X,Y,Z):
    if len(X) != len(Y) or len(X) != len(Z):
        print("Length mismatch")
        return -1
    numNodes = len(X)
    fout_vtk.write(str(numNodes)+" 3\n")
    for i in range(numNodes):
        fout_vtk.write(str(X[i])+" "+str(Y[i])+" "+str(Z[i])+"\n")


def writeConnToVTK(fid,conn):
    numConn = len(conn)
    #write header for conn
    fid.write("POLYGONS "+str(numConn)+" "+str(4*numConn)+"\n")
    for ele in conn:
        ele = list(ele) #Set does not allow indexing
        fid.write("3 "+str(ele[0])+" "+str(ele[1])+" "+str(ele[2])+"\n")

def writeConnToConnDat(fid,conn):
    numConn = len(conn)
    #write header for conn
    fid.write(str(numConn)+" Triangles\n")
    for ele in conn:
        ele = list(ele) #Set does not allow indexing
        fid.write("3 "+str(ele[0])+" "+str(ele[1])+" "+str(ele[2])+"\n")

def boundaryNodes(filename):
    [X,Y,Z] = loadNodesFromVTK(filename)
    Zmin = min(Z)
    Zmax = max(Z)
    Bndry1 = []
    Bndry2 = []
    for i in range(len(X)):
        if abs(Z[i]-Zmin) <1e-6:
            Bndry1.append(i)
        elif abs(Z[i]-Zmax) < 1e-6:
            Bndry2.append(i)
    return [Bndry1, Bndry2]

def findInteriorNodes(conn, bndry):
    """ The function takes in a connectivity table 
    and list of boundary nodes """
    interior = set([])
    for ele in conn:
        for idx in ele:
            if idx in bndry:
                setEle = set(ele)
                # indices not 'idx' (these are interior nodes)
                otherIndices = setEle - {idx}
                # consider only indices not in boundary
                innerIndices = otherIndices - (otherIndices & set(bndry)) # & = intersection 
                interior = interior | innerIndices # union
    return list(interior)

def pairUp(ghost, nodes, X, Y, Z):
    pairs = []
    for i in ghost:
        for j in nodes:
            if abs(X[i]-X[j]) <1e-4 and abs(Y[i]-Y[j])<1e-4:
                pairs.append([i,j])
    return pairs

def writeGhostInteriorPairs(filename, D):
    fid = open(filename,'w')
    # Loop over the contents of the dictionary D
    for d in D:
        fid.write("%d %d\n" %(d[0], d[1]))
    fid.close()

    
def createGhostNodePosition(idx, conn, X, Y, Z, bndry, D):
    for ele in conn:
        # If ele has an interior i1 index and a boundary b1 index
        if set(ele) & {idx} != set() and set(ele) & set(bndry) != set():
            x = list(set(ele) & set(bndry))[0] # A boundary node
            DeltaX_ix = [X[idx]-X[x], Y[idx]-Y[x], Z[idx]-Z[x]]
            Xnew = [DeltaX_ix[0]+X[D[x]], \
                    DeltaX_ix[1]+Y[D[x]], \
                    DeltaX_ix[2]+Z[D[x]]]
            return Xnew

#---------------------------------------- Script begins ------------------------------
[X,Y,Z]=loadNodesFromVTK(fileCyl)
conn = loadConnFromVTK(fileCyl)
fid_out = open(fileBump,"w")
#compute cylindrical coordinates
for i in range(len(X)):
    theta = atan2(Y[i],X[i]);
    z = Z[i];
    rho = eps*exp(-(theta**2+z**2)/(2*.4**2))
    X[i] = X[i]*(1+rho)
    Y[i] = Y[i]*(1+ rho)

writeXYZToVTK(fid_out,X,Y,Z)
writeConnToVTK(fid_out,conn)
fid_out.close()

[X,Y,Z] = loadNodesFromVTK(fileBump)
conn=loadConnFromVTK(fileBump)
[b1, b2] = boundaryNodes(fileBump)
D = pairUp(b1,b2, X, Y, Z)
Dprime = [[d[1],d[0]] for d in D]
numNodes = len(X)
i1 = findInteriorNodes(conn, b1)
i2 = findInteriorNodes(conn, b2)
D1 = [ [i1[i], numNodes+i] for i in range(len(i1))]
D1 = D1
D2 = [ [i2[i], numNodes+len(i2)+i] for i in range(len(i2))]
D2 = D2
# Two dictionaries to look up connectivites
Da = dict(D+D1)
Db = dict(Dprime+D2)
newconn = []
for ele in conn:
    # if ele has a boundary b1 node
    if (set(ele) & set(b1)) != set():
        newconn.append([Da[ele[0]], Da[ele[1]], Da[ele[2]]])
    # if ele has a boundary b2 node
    if (set(ele) & set(b2)) != set():
        newconn.append([Db[ele[0]], Db[ele[1]], Db[ele[2]]])

Xnew = []; Ynew=[]; Znew=[];
for idx in i1:
    xyz = createGhostNodePosition(idx, conn, X,Y,Z, b1, Da)
    Xnew.append(xyz[0])
    Ynew.append(xyz[1])
    Znew.append(xyz[2])

for idx in i2:
    xyz = createGhostNodePosition(idx, conn, X,Y,Z, b2, Db)
    Xnew.append(xyz[0])
    Ynew.append(xyz[1])
    Znew.append(xyz[2])

fid_out = open(fileBumpGhost,"w")
fid_nodes = open(filenodesdat,"w")
fid_conn = open(fileconndat,"w")
writeXYZToVTK(fid_out, X+Xnew, Y+Ynew, Z+Znew)
writeConnToVTK(fid_out, conn+newconn)
writeXYZToNodesDat(fid_nodes, X+Xnew, Y+Ynew, Z+Znew)
writeConnToConnDat(fid_conn, conn+newconn)
fid_out.close()
fid = open(prefix+"ghostNodeInfo.txt", 'w')
fid.write("Number of original nodes = %d \n" %len(X))
fid.write("Number of ghost nodes = %d \n" %len(Xnew))
fid.write("Ghost node indicies are %d to %d \n" %(len(X), len(X)+len(Xnew)-1))
fid.write("----------------------------\n")
fid.write("Number of original elements = %d \n" %len(conn))
fid.write("Number of ghost elements = %d \n" %len(newconn))
writeGhostInteriorPairs(prefix+"ghostBoundary1Pair.txt", D1)
writeGhostInteriorPairs(prefix+"ghostBoundary2Pair.txt", D2)
writeGhostInteriorPairs(prefix+"BoundaryPair.txt", D)
fid.close()
fid_nodes.close()
fid_conn.close()
