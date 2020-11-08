# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 17:44:01 2020

@author: Vahid Haseltalab
"""

import math, cv2
from math import *
import numpy as np
import shapely.geometry as sg
import shapely.affinity as sa
import shapely.ops as so
from PIL import Image, ImageDraw, ImageFont
from stl import mesh
from scipy import spatial
from scipy import interpolate



### Generating polygon
def polygon(r, N, angle = 0):
    angle = np.pi*angle/180
    circle = r * np.exp(1j * (np.linspace(0, 2 * math.pi, N + 1) + angle))
    out = [[round(a.real,6),round(a.imag,6)] for a in circle]
    if out[0] != out[-1]: out[-1] = out[0]
    return out
### Scale the polygon  
def scale(x,y,sc):
    return x*sc,y*sc
### Scale the polygon toward x direction  
def scalex(x,y,sc):
    x1 = [x[i]*sc for i in range(len(x))]
#    y1 = [y[i]*sc for i in range(len(y))]
    return x1,y
### Locate the polygon into center of the screen
def center(x,y,width,length):
    B = np.c_[x,y]
    if len(B) < 3: poly = sg.LineString(B)
    else: poly = sg.Polygon(B)
    minx,miny,maxx,maxy = poly.bounds
    #dis = ((width/2)-(maxx+minx)/2,(length/2)-(maxy+miny)/2)
    dis = width/2,length/2
    x = (np.array(x) + dis[0]).tolist()
    y = (np.array(y) + dis[1]).tolist()
    return x,y
### polygon operation
def PolyOp(A,B,typ = '+'):
    if typ == '+':
        poly = sg.Polygon(A).union(sg.Polygon(B))
    else:
        poly = sg.Polygon(A).difference(sg.Polygon(B))
    return list(poly.exterior.coords)
### convert polygon to binary.image
def Poly2Img(A):
    A = np.array(A)
    ### Get new coordinates
    def NewCoord(x,y,width,height):
        x,y = scale(x,y,11)
        x,y = center(x,y,width,height)
        return x,y
    ### create bitmap images
    x = A[:,0]; y = A[:,1]
    width = 1000
    height = 1000
    x,y = NewCoord(x,y,width,height)
    poly1 = sg.Polygon([[a,b] for a,b in zip(x,y)]).buffer(0)
    array = np.zeros([height, width, 3], dtype=np.uint8)
    array[:,:] = [0, 0, 0] ### defines a black pixel
    minx,miny,maxx,maxy = poly1.bounds
    for i in range(int(miny)-1,int(maxy)+1):
        l = sg.LineString([(int(minx-3),i),(int(maxx+3),i)])
        if l.intersects(poly1):
            intrsct = l.intersection(poly1)
            if type(intrsct) == sg.MultiLineString:
                intrsct = so.linemerge(intrsct)
                if type(intrsct) == sg.LineString:
                    a = list(intrsct.coords)
                    if len(a) == 3: del(a[1])
                else:
                    a = []
                    for bb in intrsct: a.extend(bb.coords)
            elif type(intrsct) == sg.LineString: a = list(intrsct.coords)
            elif type(intrsct) == sg.Point: a = []
            else: 
                intrsct = [n for n in intrsct if type(n) != sg.Point]
                a = []
                for bb in intrsct: a.extend(bb.coords)
            if len(a) == 1:
                array[i,int(a[0][0])].fill(255)
            elif len(a) == 2: array[i,int(a[0][0]):int(a[1][0])].fill(255)
            else:
                for j in range(len(a)):
                    if j%2==0:
                        array[i,int(a[j][0]):int(a[j+1][0])].fill(255)
    img = Image.fromarray(array)
    return img,array
### convert binary.image to polygon
def Img2Poly(M):
    if type(M) == str: im = cv2.imread('%s'%M)
    else: im = M
    imgray = cv2.cvtColor(im,cv2.COLOR_BGR2GRAY)
    ret,thresh = cv2.threshold(imgray,127,255,0)
    contours, hierarchy = cv2.findContours(thresh,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    coords = []
    for a in contours:
        epsilon = 0.001*cv2.arcLength(a,True)
        approx = cv2.approxPolyDP(a,epsilon,True)
        coords.append([tuple(b[0]) for b in approx])
    def percent(a,b):
        return (b.area/a.area*100)-100
    newcoords = []
    for a in coords:
        r = 20; p = 10
        if len(a) < 3: newcoords.append(a)
        else:
            while p > 1 and r > 0:
                polyg1 = sg.Polygon(a)
                polyg = sg.Polygon(a).buffer(r,join_style = 1).buffer(-r,join_style = 1)
                p = percent(polyg1,polyg)
                r -= 2
            newcoords.append(list(polyg.exterior.coords))
    return newcoords
### Simplify the polygon using a tolerance
def MinPoly(A,tol):
    poly = sg.Polygon(A)
    return list(poly.simplify(tol).exterior.coords)

### split a polygon by curve 
def PolyDisect(A,line):
    poly = sg.Polygon(A)
    line = sg.LineString(line)
    if poly.intersects(line):
        segments = so.split(poly,line)
        if segments[0].centroid.coords.xy[0][0]<segments[1].centroid.coords.xy[0][0]:
            s1 = list(segments[0].exterior.coords); s2 = list(segments[1].exterior.coords)
        else: s1 = list(segments[1].exterior.coords); s2 = list(segments[0].exterior.coords)
        return s1,s2
    else:
        print("No intersections!")
        return A
### generate polygon offset from a list of polygons, number of offsets and distance
def Offset(A,d):
    X = []; Y = []; intr = []
    if type(A) == sg.MultiPolygon:
        if A.area < 1: return [],[],None
    if type(A) != sg.MultiPolygon:
        a = sg.Polygon(A)
    else: a = A
    bp = a.buffer(d)
    if bp.area == 0: return [],[],bp
    if type(bp) == sg.multipolygon.MultiPolygon:
        for c in bp:
            c = sg.Polygon(MinPoly(c,0.05))
            x1,y1 = list(c.exterior.coords.xy)
            X.append(list(x1)); Y.append(list(y1))
            interior = []
            interi = c.interiors
            for inter in interi:
                interior.append(inter.coords.xy)
            if interior:
                intr.append(interior[0])
    else:
        if type(bp.exterior) == sg.polygon.LinearRing:
            bp = sg.Polygon(MinPoly(bp,0.05))
            x1,y1 = list(bp.exterior.coords.xy)
            X.append(list(x1)); Y.append(list(y1))
            interior = []
            interi = bp.interiors
            for inter in interi:
                interior.append(inter.coords.xy)
            if interior:
                intr.append(interior[0])
    for aa in intr:
        X.append(aa[0]); Y.append(aa[1])
    return X,Y,bp

### Slicing an STL file format. It requires STL address, and the height of slice plane
class Slice:
    def vertices(self,a):
        # importing ths stl file
        msh = mesh.Mesh.from_file('%s'%a)
        zvalues = [] ### stores the z value of a vertex 
        ### categorizing the vertices into a list based on their faces
        vrt = [[] for s in range(len(msh))] 
        for i in range(len(msh)):
            p1 = (msh[i][:3]).tolist()
            p2 = (msh[i][3:6]).tolist()
            p3 = (msh[i][6:]).tolist()
            zvalues.append(p1[-1]);zvalues.append(p2[-1])
            zvalues.append(p3[-1])
            vrt[i] = [p1,p2,p3]
        return vrt,zvalues
    
    def __init__(self,addr = 'bunny.stl',direction = 'xy'):
        if type(addr) == str:
            self.vrt,self.zvalues = self.vertices(addr)
            self.direction = direction
        else:
            self.direction = direction
            self.vrt = addr
            if direction == 'xy':
                self.zvalues = addr.reshape([len(addr)*3,3])[:,2]
            elif direction == 'xz':
                self.zvalues = addr.reshape([len(addr)*3,3])[:,1]
            elif direction == 'yz':
                self.zvalues = addr.reshape([len(addr)*3,3])[:,0]
     ### finds intersection coordinates between a point and a line 
    def eqn1(self,p1,p2,z):
        if p1[2]==p2[2]: return tuple([p1[0],p1[1]])
        ### the ratio that can apply to distance of x and y 
        t = (z-p1[2]) / (p2[2]-p1[2]) 
        return (p1[0] + (p2[0]-p1[0])*t , p1[1] + (p2[1]-p1[1])*t)
    def eqn2(self,p1,p2,z):
        if p1[1]==p2[1]: return tuple([p1[0],p1[2]])
        ### the ratio that can apply to distance of x and y 
        t = (z-p1[1]) / (p2[1]-p1[1]) 
        return (p1[0] + (p2[0]-p1[0])*t , p1[2] + (p2[2]-p1[2])*t)
    def eqn3(self,p1,p2,z):
        if p1[0]==p2[0]: return tuple([p1[1],p1[2]])
        ### the ratio that can apply to distance of x and y 
        t = (z-p1[0]) / (p2[0]-p1[0]) 
        return (p1[1] + (p2[1]-p1[1])*t , p1[2] + (p2[2]-p1[2])*t)
    ### checks whether the z plane is crossing through the line
    def checkline(self,zl,z):
        if z <= np.max(zl) and z >= np.min(zl):
            return True
        else: return False
    ### finds intersection coordinates between a plane and a triangular facet        
    def trintersct(self,l,z):
        l = np.array(l)
        if self.direction == 'xy' and (l[:,2] == z).all(): return []
        elif self.direction == 'xz' and (l[:,1] == z).all(): return []
        elif self.direction == 'yz' and (l[:,0] == z).all(): return []
        inlst = []
        for i in range(3):
            pt1 = l[i];pt2 = l[i-1]
            if self.direction == 'xy':
                zl = [pt1[2],pt2[2]]
                if self.checkline(zl,z):
                    p = self.eqn1(pt1,pt2,z)
                    inlst.append(p)
            elif self.direction == 'xz':
                zl = [pt1[1],pt2[1]]
                if self.checkline(zl,z):
                    p = self.eqn2(pt1,pt2,z)
                    inlst.append(p)
            elif self.direction == 'yz':
                zl = [pt1[0],pt2[0]]
                if self.checkline(zl,z):
                    p = self.eqn3(pt1,pt2,z)
                    inlst.append(p)
        if len(inlst) == 3:
            if np.array_equal(np.round(inlst[0],4),np.round(inlst[1],4)): del(inlst[1])
            else: del(inlst[2])
        if inlst:
            if np.array_equal(np.round(inlst[0],4),np.round(inlst[1],4)): return []
        return inlst
    #The function clear will remove empty lists from S
    def clear(self,A):
        n=len(A);s=0;k=0
        for i in range(n):
            s+=int(not A[i])
        W=[[0 for i in range(1)]for j in range(n-s)]
        for i in range(n):
            if A[i]:
                W[k]=A[i]
                k+=1
        return W
    def totuple(self,a): ## convert an array to nested tuple
        try:
            return tuple(self.totuple(i) for i in a)
        except TypeError:
            return a
    def order(self,pnts,prec = 5):
        if not pnts: return []
        pnts = np.round(pnts,prec)
        pnts = self.totuple(pnts)
        wires = [];Q1 = set()
        for pn in pnts:
            #pn = simplify(pn,0.01)
            if pn not in Q1:
                Q = list(pn)
                Q1.add(pn)
                a1 = Q[0]; a2 = Q[1]; a3 = Q[-2]; a4 = Q[-1]
                c = 0
                while c < 1:
                    run = False
                    for b in pnts:
                        #b = simplify(b,0.01)
                        if b not in Q1:
                            b1 = b[0]; b2 = b[-1]
                            if a4 == b1:
                                run = True
                                Q = Q + list(b)
                                Q1.add(b)
                                a3 = b1; a4 = b2
                            elif a1 == b2:
                                run = True
                                Q = list(b) + Q
                                Q1.add(b)
                                a1 = b1; a2 = b2
                            elif a4 == b2:
                                run = True
                                Q = Q + list(b)[::-1]
                                Q1.add(b)
                                a3 = b2; a4 = b1
                            elif a1 == b1:
                                run = True
                                Q = list(b)[::-1] + Q
                                Q1.add(b)
                                a1 = b2; a2 = b1
                    if not run: c+=1
                wires.append(Q)
        if len(wires) == 1:
            return wires[0]
        if prec == 0: return wires[0]
        else: return self.order(wires,prec-1)
    def cross_section(self,z):
        z1 = np.array(self.zvalues).reshape((int(len(self.zvalues)/3),3))
        idx = np.where((z1>z-5)&(z1<z+5))[0]
        fnum = set()
        L = []
        for a in list(set(idx)):
            if a not in fnum:
                fnum.add(a)
                cr = self.trintersct(self.vrt[a],z)
                if cr: L.append(cr)
        L1 = np.round(L,3)
        lines = list(so.polygonize(L1.tolist()))
        if lines:
            return list(lines[0].exterior.coords)
        else:
            return list(self.order(L))
### combine the subparts in 2D
def combine(plys,boolean_list):
    if plys:
        base = 0
        for a in plys:
            if a.area!= 0:
                base = a
                break
        if not base.is_valid and base != 0: base = base.buffer(0)
        for j in range(1,len(plys)):
            if plys[j].area > 0.5:
                if boolean_list[j] == 1:
                    base = base.union(plys[j])
                elif boolean_list[j] == -1:
                    base = base.difference(plys[j])
        try:
            bs = [a for a in base if a.area>2]
        except: bs = [base]
        return bs
    else: return []
### Gives cross_section of combines subparts
def xsect(parts,b_list,z):
    slices1 = [Slice(parts[i],'xy') for i in range(len(parts))]
    sec1 = [b.cross_section(z) for b in slices1]
    plys1 = [sg.Polygon(a).buffer(0) for a in sec1]
    base1 = combine(plys1,b_list)
    lst1 = [list(a.exterior.coords) for a in base1]
    return lst1

### 2D transformation
def HTM2D(A,alpha,x0,y0):
    poly = sg.Polygon(A)
    poly = sa.rotate(poly,alpha)
    poly = sa.translate(poly,x0,y0)
    return list(poly.exterior.coords)
### apply 3D transformation on an object
def Transform(obj,thx=0,thy=0,thz=0,dx=0,dy=0,dz=0):
    thx = radians(thx); thy = radians(thy); thz = radians(thz)
    col1 = (cos(thy)*cos(thz),-cos(thy)*sin(thz),sin(thy),dx)
    col2 = (sin(thx)*sin(thy)*cos(thz)+cos(thx)*sin(thz),-sin(thx)*sin(thy)*sin(thz)+cos(thx)*cos(thz),-sin(thx)*cos(thy),dy)
    col3 = (-cos(thx)*sin(thy)*cos(thz)+sin(thx)*sin(thz),cos(thx)*sin(thy)*sin(thz)+sin(thx)*cos(thz),cos(thx)*cos(thy),dz)
    col4 = (0,0,0,1)
    T = np.column_stack([col1,col2,col3,col4])
    obj1 = np.insert(obj,3,1,axis=1)
    result = np.matmul(obj1,T)
    result = np.delete(result,3,axis=1)
    return np.round(result,5)

### generating 2D triangular mesh (random points are distributed based on the value of percentage from 0 to 1)
def mesh2D(poly):
    pp = sg.Polygon(poly)
    xmin,ymin,xmax,ymax = pp.bounds
    points = np.array(poly)
    if xmax-xmin > ymax-ymin: dis = (xmax-xmin)/5
    else: dis = (ymax-ymin)/5
    x,y = np.meshgrid(np.arange(xmin-1, xmax+1, dis/8), np.arange(ymin-1, ymax+1, dis/8))
    pgrid = np.c_[x.ravel(),y.ravel()]
    for p in pgrid:
        p11 = sg.Point(p)
        if p11.within(pp) and min(np.linalg.norm(points-p,axis = 1)) > dis:
            points = np.vstack([points,p])
    ppoints = sg.MultiPoint(points)
    tri = so.triangulate(ppoints)
    ff = [a for a in tri if a.within(pp)]
    vert = [list(b.exterior.coords) for b in ff]
    return vert, points
### increase the number of points on a 2D curve
def ssample(poly,p):
    X = np.array([a[0] for a in poly]); Y = np.array([a[1] for a in poly])
    s=np.insert(np.cumsum(np.absolute(np.add(np.diff(X),1j*np.diff(Y)))),0,0)
    if p<1:p=np.ceil(s[-1]/p)
    n=len(X);delta=s[-1]/p
    xs=np.zeros((p+1,1));ys=np.zeros((p+1,1))
    xs[0]=X[0];ys[0]=Y[0]
    for i in range(1,p+1):
        sr=i*delta;j=np.sum(s<sr)
        if j!=n:
            u=(sr-s[j-1])/(s[j]-s[j-1])
            dx=X[j]-X[j-1];dy=Y[j]-Y[j-1]
            xs[i]=X[j-1]+u*dx;ys[i]=Y[j-1]+u*dy
        else: xs[i]=X[-1];ys[i]=Y[-1]
    return np.asarray([[a[0],b[0]] for a,b in zip(xs,ys)])
### The least common multiple (L.C.M.) of two numbers
def lcm(x, y):
   if x > y:
       greater = x
   else:
       greater = y
   while(True):
       if((greater % x == 0) and (greater % y == 0)):
           lcm = greater
           break
       greater += 1
   return lcm
### prepare scale list
def scale_(scale,n):
#    length = np.cumsum(np.linalg.norm(np.diff(path,axis=0),axis=1))[-1]
    dis = int(n/(len(scale)-1))
    scl = np.ones(n)
    for i in range(len(scale)-1):
        scl[i*dis:(i+1)*dis] = np.linspace(scale[i],scale[i+1],dis).tolist()
    scl[-(n%dis):] = scl[-(n%dis)-1]
    return scl
### polygon morphing, n is the number of polygons in between
### twist angle (n) and scaling factors will also affect the morphed polygons
def PlygnMorph(poly1,poly2,n,angle,scl):
    scl = scale_(scl,n)
    lst = range(n)
    if len(poly1) != len(poly2):
        num = lcm(len(poly2)-1, len(poly1)-1)
        if num < 8: num *= 2
        if num > 25:
            num = max([len(poly1),len(poly2)])*2
        polygn1 = ssample(poly1,num); polygn2 = ssample(poly2,num)
    else:
        polygn1 = poly1; polygn2 = poly2
    def midply(poly1,poly2):
        newpoly = []
        for i in range(len(poly1)):
            if list(poly1[i]) == list(poly2[i]): newpoly.append(poly1[i])
            else:
                ln = sg.LineString([poly1[i],poly2[i]])
                newpoly.extend(list(ln.interpolate(ln.length/2).coords))
        return newpoly
    polys = [np.array(polygn1).tolist(),np.array(polygn2).tolist()]
    seq = [0,lst[-1]]
    running = True
    while running:
        running = False
        Q = [a for a in seq]; pl = [b for b in polys]
        for i in range(len(Q)-1):
            m = lst[Q[i+1]] - lst[Q[i]]
            if m > 1:
                running = True
                polys.insert(polys.index(pl[i])+1,midply(pl[i],pl[i+1]))
                seq.insert(seq.index(Q[i])+1,int((Q[i]+Q[i+1])/2))
    plys = [sg.Polygon(a) for a in polys]
    result = []
    for i in range(len(plys)):
        cen = np.round(list(plys[i].centroid.coords),4)[0]
        p = sa.scale(plys[i],xfact=scl[i],yfact=scl[i])
        p = sa.rotate(p,angle/(n-1)*i,sg.Point(cen))
        result.append(list(p.exterior.coords))
    return np.array(result)
### transformation matrix
def Matrix(v1,p):
    n = np.linalg.norm(v1)
    a = [0,0,1]
    b = v1/n
    v = np.cross(a,b)
    c = np.dot(a,b)
    s = np.linalg.norm(v)
    if not np.any(v): s = 0.00000000000001
    I = np.identity(3)
    vXStr = '{} {} {}; {} {} {}; {} {} {}'.format(0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0)
    k = np.matrix(vXStr)
    r = I + k + np.matmul(k,k) * ((1 -c)/(s**2))
    r = np.insert(r,3,p,axis=1)
    r = np.insert(r,3,[0,0,0,1],axis=0)
    return r

### find the rotation angles between two vectors
def angles2(v):
    n = np.linalg.norm(v); un = v/n
    thx = 0
    thy = np.arccos(un[2])
    thz = np.arctan2(un[1],un[0])
    return np.degrees(thx),np.degrees(thy),np.degrees(thz)
### map 2D polygons to 3D using a path
def mapto3D11(polys,path,typ='normal'):
    plys = []
    thx,thy,thz = angles2(np.array(path[1])-np.array(path[0]))
    dx,dy,dz = np.array(path[0])
    plys.append(Transform(polys[0],thx,thy,thz,dx,dy,dz))
    #agl = []
    for i in range(len(path)-1):
        v = np.array(path[i])-np.array(path[i-1])
        thx,thy,thz = angles2(v)
        if (path[0] == path[-1]).all() and i == len(path)-2:
            thx,thy,thz = 0,0,0
        if typ =='normal': dx,dy,dz = np.array(path[i+1])
        else: dx,dy,dz = 0,0,0
        b = Transform(polys[i],thx,thy,thz,dx,dy,dz)
        if isinstance(b,list): b = np.array(b)
        plys.append(b)
    return plys
### map 2D polygons to 3D using a path
def mapto3D(polys,path,typ='normal'):
    plys = []
    vec, path = redefine_path(path)
    if typ == 'normal': T = Matrix([0,0,1],path[0])
    else: T = Matrix([0,0,1],[0,0,0])
    if typ == 'normal' and (path[0] == path[-1]).all():
        T = Matrix([0,0,1],[0,0,0])
    obj1 = np.insert(polys[0],3,1,axis=1)
    result = np.array(np.matmul(T,obj1.transpose()))
    plys.append(np.c_[result[0],result[1],result[2]])
    c = 0
    for i in range(1,len(path)):
        c+=1
        if typ =='normal': T = Matrix(vec[i-1],path[i])
        else: T = Matrix(vec[i-1],[0,0,0])
        if (path[i] == path[i-1]).all() and i != 0: c-=1
        obj1 = np.insert(polys[c],3,1,axis=1)
        result = np.array(np.matmul(T,obj1.transpose()))
        plys.append(np.c_[result[0],result[1],result[2]])
    if typ != 'normal': T = Matrix([0,0,1],[0,0,0])
    else: T = Matrix(vec[-1],path[-1])
    if typ == 'normal' and (path[0] == path[-1]).all():
        T = Matrix([0,0,1],[0,0,0])
    obj1 = np.insert(polys[-1],3,1,axis=1)
    result = np.array(np.matmul(T,obj1.transpose()))
    plys.append(np.c_[result[0],result[1],result[2]])
    return plys, path
### interpolate 3D curve and increase its number of points (p); line is a numpy array
def ssample3D(line,p):
    line = np.array(line)
    if (line==line[0]).all():
        return np.array([line[0]]*(p+1))
    X = line[:,0]; Y = line[:,1];Z = line[:,2]
    s=np.insert(np.cumsum(np.linalg.norm(np.diff(line,axis=0),axis=1)),0,0)
    if p<1:p=np.ceil(s[-1]/p)
    n=len(X);delta=s[-1]/p
    xs=np.zeros((p+1,1));ys=np.zeros((p+1,1));zs = np.zeros((p+1,1))
    xs[0]=X[0];ys[0]=Y[0]; zs[0]=Z[0]
    for i in range(1,p+1):
        sr=i*delta;j=np.sum(s<sr)
        if j!=n:
            u=(sr-s[j-1])/(s[j]-s[j-1])
            dx=X[j]-X[j-1];dy=Y[j]-Y[j-1];dz = Z[j]-Z[j-1]
            xs[i]=X[j-1]+u*dx;ys[i]=Y[j-1]+u*dy;zs[i]=Z[j-1]+u*dz
        else: xs[i]=X[-1];ys[i]=Y[-1];zs[i]=Z[-1]
    return np.c_[xs,ys,zs]
### find indices of 2D mesh with respect to its polygon
def nodes(poly,mesh):
    nds = []
    if isinstance(poly[0],list) or isinstance(poly,np.ndarray): poly = [tuple(a) for a in poly]
    for a in mesh:
        q = []
        for b in a:
            q.append(poly.index(b))
        nds.append(q)
    return nds
### save vertices and faces as stl file with given name
def savestl(name,vertices,faces):
    f = open('%s.stl'%name,'w')
    f.write('solid %s\n'%name)
    for a in vertices[faces]:
        ln1,ln2 = np.diff(a,axis=0)
        cross = np.cross(ln1,ln2)
        normal = cross/np.linalg.norm(cross)
        f.write('\tfacet normal %f %f %f\n'%(normal[0],normal[1],normal[2]))
        f.write('\t\touter loop\n')
        for j in range(3):
            f.write('\t\t\tvertex %f %f %f\n'%(a[j][0],a[j][1],a[j][2]))
        f.write('\t\tendloop\n')
        f.write('\tendfacet\n')
    f.write('endsolid %s'%name)
    f.close()
### create obj file 
def saveobj(name,vertices,faces):
    f = open('%s.obj'%name,'w')
    for a in vertices:
        f.write('v %s %s %s\n'%(a[0],a[1],a[2]))
    for a in faces:
        f.write('f %s %s %s\n'%(a[0]+1,a[1]+1,a[2]+1))
    f.close()
### sweep function which takes two profiles and a guide curve
### as well as twist angle and scaling factors (as a list)
def sweep(poly1,poly2,path,TwistAngle=0,scale1=None):
    path = np.array(path)
    resolution = 3
    if scale1 == None: scale1 = [1,1]
    poly1 = np.round(np.array(poly1),5); poly2 = np.round(np.array(poly2),5)
    if poly1.tolist() == poly2.tolist() and (path[0] == path[-1]).all():
        if np.count_nonzero(np.count_nonzero(poly1,axis=1) == 0) > 2:
            poly = opt_line(poly1[:-1],resolution)
        else: poly = opt_line(poly1,resolution)
        lsts = np.array([poly]*len(path))
        lsts1 = np.insert(np.array(lsts), 2, 0, axis=2)
        curves2,path1 = mapto3D(lsts1,path,'rot'); path2 = path1
    else:
        path1 = interplt(path,resolution,int(len(path)*4))
        if poly1.tolist() == poly2.tolist() and TwistAngle==0.0 and scale1==scale1[0]: lsts = np.array([poly1]*len(path1))
        else: lsts = np.array(PlygnMorph(poly1,poly2,len(path1),TwistAngle,scale1))
        lsts1 = np.insert(lsts, 2, 0, axis=2)
        polys,path2 = mapto3D(lsts1,path1)
        sch1 = np.max(np.linalg.norm(poly1-np.array([0,0]),axis=1))
        sch2 = np.max(np.linalg.norm(poly2-np.array([0,0]),axis=1))
        searchradius = np.max([sch1,sch2])+5
        tree = spatial.cKDTree(path1)
        curves = [[] for _ in range(len(polys[0]))]
        for j in range(len(polys[0])): curves[j].append(polys[0][j])
        for i in range(1,len(path2)):
            idx = []; pidx = tree.query_ball_point(path2[i],searchradius)
            a = [np.linalg.norm(polys[i]-path1[k],axis=1) for k in pidx]
            rf = a.pop(pidx.index(np.where((path1 == path2[i]).all(axis=1))[0][0]))
            idx = np.where(np.greater(np.array(a),rf).all(axis=0))[0]
            for j in idx: curves[j].append(polys[i][j])
        smoothness = 0.5 + TwistAngle * 0.003 
        size = len(curves[0])*2
        for a in curves:
            if len(a) > size: size = len(a)
        curves2 = [interplt(a,resolution,size,smoothness) for a in curves]
    num_curve = len(curves2[1])
    vertices = np.concatenate(np.array(curves2))
    faces = []
    for i in range(len(curves2)-1):
        for j in range(num_curve-1):
            id1 = (num_curve)*i+j
            id2 = (num_curve)*i+j+1
            id3 = (num_curve)*(i+1)+j
            id4 = (num_curve)*(i+1)+j+1
            faces.append([id3,id2,id1]); faces.append([id3,id4,id2])
    if (path1[0] == path1[-1]).all():
        fc = faces
    else: ## closing the ends
        mesh1,poly1 = mesh2D(lsts[0])
        mesh2, poly2 = mesh2D(lsts[-1])
        newplys = lsts1[1:-1].tolist()
        newplys = [np.insert(poly1, 2, 0, axis=1).tolist()]+newplys+[np.insert(poly2, 2, 0, axis=1).tolist()]
        mesh1 = [a[::-1] for a in mesh1]
        newpolys,_ = mapto3D(newplys,path2)
        firstf = nodes(poly1,mesh1); firstf = np.delete(np.array(firstf),3,1)
        secondf = nodes(poly2,mesh2); secondf = np.delete(np.array(secondf),3,1)
        num_vert = len(vertices)
        vertices = np.append(vertices,np.array(newpolys[-1]),0)
        vertices = np.append(vertices,np.array(newpolys[0]),0)
        secondf = secondf+num_vert
        firstf = firstf+(num_vert+len(poly2))
        fc = np.array(firstf.tolist()+faces+secondf.tolist())
    return vertices,fc
### interpolate a 3D curve with given number of points (p)
def interplt(curve,res,p,smoothness=1.5):
    vec = np.diff(curve,axis=0); gd = np.array(grd(vec))
    tol = 2  ## bounds the end points with different interpolation curve
    if curve[0].tolist() == curve[-1].tolist(): return ssample3D(curve,p+20)
    if len(curve) < 10: return opt_line(curve,res)
    curve = np.array(curve)
    _, idx = np.unique(curve,axis=0,return_index=True)
    curve = curve[np.sort(idx)]
    x = curve[:,0]; y = curve[:,1]; z = curve[:,2]
    if min(gd) > 0.99:
        tck, u = interpolate.splprep([x,y,z], k=3, s=0)
        u_fine = np.linspace(0,1,p+20)
        x1,y1,z1 = interpolate.splev(u_fine, tck)
        return np.c_[x1,y1,z1]
    if (not x.any() and not y.any()) or (not x.any() and not z.any()) or (not y.any() and not z.any()):
        return opt_line(curve,res)
    l, r = [(1, (0, 0, 0))], [(2, (0, 0, 0))]
    tck, u = interpolate.splprep([x[tol:-tol],y[tol:-tol],z[tol:-tol]], k=3, s=smoothness, t = curve[tol:-tol])
    u_fine = np.linspace(0,1,p)
    x1,y1,z1 = interpolate.splev(u_fine, tck)
    x_old0 = np.concatenate([x[:tol],x1[:tol]],axis=0); y_old0 = np.concatenate([y[:tol],y1[:tol]],axis=0);z_old0 = np.concatenate([z[:tol],z1[:tol]],axis=0)
    x_old1 = np.concatenate([x1[-tol:],x[-tol:]]); y_old1 = np.concatenate([y1[-tol:],y[-tol:]]); z_old1 = np.concatenate([z1[-tol:],z[-tol:]])
    tckp, u = interpolate.splprep([x_old0,y_old0,z_old0], s=3, k=3)
    clamped_spline = interpolate.make_interp_spline(u, np.array([x_old0, y_old0, z_old0]).T, bc_type=(l, r))
    xnew0, ynew0, znew0 = clamped_spline(np.linspace(0,1,len(x_old0)*3)).T
    tckp, u = interpolate.splprep([x_old1,y_old1,z_old1], s=3, k=3)
    clamped_spline = interpolate.make_interp_spline(u, np.array([x_old1, y_old1, z_old1]).T, bc_type=(l, r))
    xnew1, ynew1, znew1 = clamped_spline(np.linspace(0,1,len(x_old1)*3)).T
    tck0, u0 = interpolate.splprep([x_old0,y_old0,z_old0], k=3, s=None)
    u_fine0 = np.linspace(0,1,len(x_old0)*3)
    xnew0,ynew0,znew0 = interpolate.splev(u_fine0, tck0)
    tck1, u1 = interpolate.splprep([x_old1,y_old1,z_old1], k=3, s=None)
    u_fine1 = np.linspace(0,1,len(x_old1)*3)
    xnew1,ynew1,znew1 = interpolate.splev(u_fine1, tck1)
    x1 = np.concatenate([xnew0,x1[tol:-tol],xnew1])
    y1 = np.concatenate([ynew0,y1[tol:-tol],ynew1])
    z1 = np.concatenate([znew0,z1[tol:-tol],znew1])
    return np.c_[x1,y1,z1]
### increase points based on the resolution
def opt_line(obj,resolution):
    obj = np.array(obj)
    if len(obj[0]) == 2:
        opt_obj = []
        for i in range(len(obj)-1):
            a = sg.LineString(obj[i:i+2]).length
            if a > resolution:
                p = int(a/resolution)+1
                opt_obj.extend(ssample(obj[i:i+2],p)[:-1])
            else:
                opt_obj.append(obj[i])
        opt_obj.append(obj[-1])
    else:
        opt_obj = []
        for i in range(len(obj)-1):
            a = np.linalg.norm(np.diff(obj[i:i+2],axis=0),axis=1)[0]
            if a > resolution:
                p = int(a/resolution)+1
                opt_obj.extend(ssample3D(obj[i:i+2],p)[:-1])
            else:
                opt_obj.append(obj[i])
        opt_obj.append(obj[-1])
    return np.array(opt_obj)

### regenerate curve and vectors based on their gradients 
def redefine_path(curve):
    curve = np.vstack((curve[0],[curve[i] for i in range(1,len(curve)) if not (curve[i] == curve[i-1]).all()]))
    vec = np.diff(curve,axis=0)
    while min(grd(vec)) < 0.9:
        for i in range(len(vec)-1):
            vec, curve = check_angle(vec,curve,i)
    return vec, curve
### add vector and point to vector/path list
def check_angle(vec,path, idx):
    v1 = vec[idx]; v2 = vec[idx+1]
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    un1 = v1/n1; un2 = v2/n2
    if np.inner(un1,un2) < 0.9:
        vec = np.insert(vec,idx+1,[(v1[0]+v2[0])/2,(v1[1]+v2[1])/2,(v1[2]+v2[2])/2], axis=0)
        path = np.insert(path,idx+1,path[idx+1],axis=0)
    return vec, path
## find gradient list of vectors
def grd(vec):
    nn = np.linalg.norm(vec,axis=1)
    unit_vec = vec/nn[:,np.newaxis]
    grdnt = []
    for i in range(len(unit_vec)-1):
        grdnt.append(np.inner(unit_vec[i],unit_vec[i+1]))
    return grdnt

### binary operations between two arrays. the typ must be choosed between "or", "and" and "not" as a string
def bin_opr(arr1,arr2,typ='or'):
    width = 1920; height = 1080
    arr3 = np.zeros([height, width, 3], dtype=np.uint8)
    arr3[:,:] = [0, 0, 0] ### defines a black pixel
    cond1 = arr1 == 255; cond2 = arr2 == 255
    if typ == 'or':
        arr3[np.where(cond1|cond2)] = 255
    elif typ == 'not':
        arr3[np.where(cond1& np.logical_not(arr2 == 255))] = 255
    elif typ == 'and':
        arr3[np.where(cond1&cond2)] = 255
    return arr3

# =============================================================================
# generate text
# =============================================================================
def text_maker(txt):
    img = Image.new('RGB', (1920, 1080), color = (0, 0, 0))
    fnt = ImageFont.truetype('cour.ttf', 200)
    d = ImageDraw.Draw(img)
    d.text((1080/2,1080/2), "%s"%txt, font=fnt, fill=(255, 255, 255))
#    img.show()
    array = np.asarray(img)
    
    ### mirror
    #img1 = img.transpose(Image.FLIP_LEFT_RIGHT)
    
    ### --------------------------------------
    
    img.save('pil_text_font.png')
    coords = Img2Poly('pil_text_font.png')
    hh = []
    for a in coords:
        a = np.array(a)
        hh.append(scale(a[:,0],a[:,1],1/8))
    gg = [list(zip(a[0],a[1])) for a in hh]
    return gg
# =============================================================================
# project text on a polygon
# =============================================================================
def project(coords,p,prj_cen,i):
    #plot_slice(gg)
    #plot_slice(list(l.coords))
    pp = sg.MultiPolygon([sg.Polygon(a) for a in coords])
    minx,miny,maxx,maxy = pp.bounds
    minx1,miny1,maxx1,maxy1 = p.bounds
    cen = (maxx+minx)/2
    pp = sa.translate(pp,-cen+prj_cen[0],-(miny-maxy1-20))
    minx,miny,maxx,maxy = pp.bounds
    if maxy > i+miny:
        l = sg.LineString([(int(minx-3),i+miny),(int(maxx+3),i+miny)])
        intrsct = l.intersection(pp)
        if type(intrsct) == sg.multilinestring.MultiLineString:
            for b in intrsct:
                points = list(b.coords)
                lines = [sg.LineString(ln) for ln in list(zip(points,[prj_cen]*2))]
                ln1,ln2 = [a.intersection(p) for a in lines]
                p1 = list(lines[0].coords)[0]; p2 = list(ln1.interpolate(1).coords)[0]
                p3 = list(ln2.interpolate(1).coords)[0]; p4 = list(lines[1].coords)[0]
                p = p.difference(sg.Polygon([p1,p2,p3,p4]))
    return p

### read coordinates from text file and convert then into np.array
def Load(filename, type1 = 'path'):
    f = open('%s.txt'%filename,'r')
    data = f.readlines()
    path = []
    for a in data:
        a = a.split(',')
        if type1 == 'path':
            path.append([float(a[1]),float(a[2]),float(a[3])])
        else:
            path.append([float(a[1]),float(a[2])])
    f.close()
    return np.array(path)

