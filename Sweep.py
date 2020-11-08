# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 17:44:01 2020

@author: Vahid Haseltalab
"""

from math import *
import numpy as np
import shapely.geometry as sg
import shapely.affinity as sa
import shapely.ops as so
from scipy import spatial
from scipy import interpolate

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
