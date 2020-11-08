# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 17:04:02 2020

@author: Vahid
"""

from LIPRO import *
boolean_list = []; parts = []

profile1 = Load('Text coordinates/turb1','profile')
path = polygon(4,50)
path = np.round(np.array(Transform(np.insert(path,2,0,axis=1),-90,0,180,0,0,0)),decimals=8)
part1 = sweep(profile1,profile1,path,0,0.5)
part11 = np.array(Transform(part1[0],90,0,0,0,0,0))
parts.append(part11[part1[1]]); boolean_list.append(1)
#savestl('prt1',part11,part1[1])
path = Load('txt/path2','path')
profile2_1 = Load('Text coordinates/profile2_1','profile'); profile2_2 = Load('Text coordinates/profile2_2','profile')
part2 = sweep(profile2_1,profile2_2,path,-50)
part22 = np.array(Transform(part2[0],0,270,0,34.8,0,10.5))
for i in range(12):
    prt = np.array(Transform(part22,0,0,i*30,0,0,0))
    parts.append(prt[part2[1]]); boolean_list.append(1)
#    savestl('prt2_%d'%i,prt,part2[1])
profile1 = polygon(4.5,100); path = np.array([[0,0,0],[0,0,40.5]])
part3 = sweep(profile1,profile1,path)
parts.append(part3[0][part3[1]]); boolean_list.append(-1)
#savestl('prt3',part3[0],part3[1])
