# LIPRO
> This is the documentation of LIPRO which is written in Python. Examples provided here are self-contained and presented for better understanding of users.

### Required libraries

[OpenCV-python](https://pypi.org/project/opencv-python "LCO")\
[Numpy](https://numpy.org/doc "LCO")\
[Shapely](https://shapely.readthedocs.io/en/stable/manual.html "LCO")\
[Numpy-stl](https://pypi.org/project/numpy-stl/ "LCO")\
[Scipy](https://www.scipy.org/docs.html "LCO")

Since LIPRO is written in Python, it uses the basic data structure list to present vectors and directions. 
```sh
v1 = [0,0,1]

curve = [[0,0,0],[10,5,5],[10,10,10]]
```
In order to have 2D polygons, function polygon(radius,Num of edges, starting angle) is presented.

```sh
poly1 = polygon(10,4,0)
```
Having a guide curve and two closed 2D polygons (start and end of the guide curve), we can generate a 3D model and present it as STL or OBJ files.

```sh
vertices, faces = sweep(poly1, poly2, guide curve, twist angle=0, scale=None)
savestl('name of the part', vertices, faces)
```

There is also a transformation function to move and rotate each part in space based on cartesian coordinates. 
```sh
part = Transform(object,thx=0,thy=0,thz=0,dx=0,dy=0,dz=0)
```
where thx,thy,thz are rotation around x,y,z axes and dx,dy,dz are translation along x,y,z axes. 


## Usage example

### Cube
To generate a cube, we only need to move a single square profile along a path as presented below. 

```python
from LIPRO import *

P1 = polygon(10,4)
G1 = [[0,0,0],[0,0,10*np.sqrt(2)]]
part1 = sweep(P1,P1,G1)

#savestl('cube',part1[0],part1[1])
```

<img src="https://user-images.githubusercontent.com/53440292/96132563-3ac15780-0f03-11eb-84f6-8301e2362cc3.png" width="30%">

### Bishop

Some profiles are too complicated to be designed by combination of polygons. Therefore, they can be loaded and used with sweep function.

```python
from LIPRO import *

profile1 = Load('text files/Bishop','profile')
path = polygon(4,50,90)
path = Transform(np.insert(path,2,0,axis=1),0,90,0,-4,0,0)
part1 = sweep(profile1,profile1,path)

#savestl('Bishop', part1[0],part1[1])
```

<img src="https://user-images.githubusercontent.com/53440292/96174975-d6b68780-0f32-11eb-953e-905a7f5765ed.JPG" width="40%">

### Donuts

Two parts can be designed separately and can join together to form a single part.

```python
from LIPRO import *

G1 = np.insert(polygon(1,50),2,0,axis=1)
G1 = Transform(G1,90,0,0,0,0,0)
P1 = np.array(polygon(1,50)) + [-1,0]
part1 = sweep(P1,P1,G1)
part1_ver = Transform(part1[0],0,0,0,0,0,1)

G2 = np.insert(polygon(1,50),2,0,axis=1)
G2 = Transform(G2,0,90,90,0,0,0)
P2 = np.array(polygon(1,50,90)) + [0,1]
part2 = sweep(P2,P2,G2)
part2_ver = Transform(part2[0],0,0,0,0,0,3)

#savestl('ring',part1_ver,part1[1])
#savestl('ring2',part2_ver,part2[1])
```

<img src="https://user-images.githubusercontent.com/53440292/96175871-2cd7fa80-0f34-11eb-917d-5444f4c313c8.jpg" width="50%">

### A path with multiple profiles

If along the guide curve there must be located more than two profiles, we can decompose the part into multiple subparts and use the sweep function to make each subpart. Afterwards, these subparts can be joined together to form the original part.
The number of subparts can be obtained from,

```sh
N = Number of profiles - 1
```
For instance in the following example there are four subparts and five profiles as, 

<img src="https://user-images.githubusercontent.com/53440292/97698472-9b32c600-1ab9-11eb-9523-9b28b2b7ce8e.png" width="50%">

```python
from LIPRO import *

d = 11
parts = []; b_list = []
G = Transform(np.insert(np.array(polygon(100,100)[:int(100/2)]),2,0,1),90,180,0,100,0,0)
G = G[:int(len(G)/4)+1]
P1 = polygon(d,3)
P2 = polygon(d,4)
part1 = sweep(P1,P2,G)
parts.append(part1[0][part1[1]]);b_list.append(1)
#savestl('m1',part1[0],part1[1])

P3 = polygon(d,5)
part2 = sweep(P2,P3,G)
part22 = Transform(part2[0],0,43.2,0,G[-1][0]-0.2,0,G[-1][2])
parts.append(part22[part2[1]]);b_list.append(1)
#savestl('m2',part22,part2[1])

P4 = polygon(d,6)
part3 = sweep(P3,P4,G)
G1 = Transform(G,0,43.2,0,G[-1][0],0,G[-1][2])
part33 = Transform(part3[0],0,86.4,0,G1[-1][0]-0.3,0,G1[-1][2])
parts.append(part33[part3[1]]);b_list.append(1)
#savestl('m3',part33,part3[1])

P5 = polygon(d,7)
part4 = sweep(P4,P5,G)
G2 = Transform(G,0,86.4,0,G1[-1][0],0,G1[-1][2])
part44 = Transform(part4[0],0,129.6,0,G2[-1][0]-0.3,0,G2[-1][2])
parts.append(part44[part4[1]]);b_list.append(1)
#savestl('m4',part44,part4[1])
```

<img src="https://user-images.githubusercontent.com/53440292/97699073-9cb0be00-1aba-11eb-95e4-fd5a30cc22af.jpeg" width="40%">
<img src="https://user-images.githubusercontent.com/53440292/97699114-aa664380-1aba-11eb-9b34-1e9278e5009d.jpeg" width="40%">

All the parts are stored into a list called **parts**, and a list of ones and mines ones defines the boolean operation of each subpart (**b_list**). Having these list of data along with the specified elevation we can join the subparts at each elevation using the following script. 

```python
cross_section = xsect(parts,b_list,90)
```
Plotting this cross section is possible using function below.

```python
def plot(lst):
    import matplotlib.pyplot as plt
    for i in range(len(lst)):
        ff = np.array(lst[i])
        plt.plot(ff[:,0],ff[:,1])
        #plt.axes().set_aspect('equal')  ## uncomment this if you need equal aspect ratio
```
The result of this plot is shown as, 

<img src="https://user-images.githubusercontent.com/53440292/97710063-4ea4b600-1acc-11eb-99f2-88eaa6278690.jpg" width="60%"> <img src="https://user-images.githubusercontent.com/53440292/97709299-2799b480-1acb-11eb-83b7-1e656033997b.png" width="30%">


Also using the following lines, multiple cross_sections of this part can be plotted in a single window.

```python
for i in range(0,102,10):
    cr = xsect(parts,b_list,i+0.1)
    plot(cr)
```

<img src="https://user-images.githubusercontent.com/53440292/97894387-7ea9be80-1d43-11eb-9c8d-7295ac5dfdf8.png" width = "80%">

