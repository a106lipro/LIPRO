# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 14:55:05 2020

@author: Vahid
"""

from LIPRO import *
boolean_list = []; parts = []
### part1
profile = Load('Text coordinates/teapot2','profile')
path = polygon(4,50)
path = np.round(np.array(Transform(np.insert(path,2,0,axis=1),-90,0,180,-4,0,0)),decimals=8)
part1 = sweep(profile,profile,path,0)
part11 = np.array(Transform(part1[0],90,0,0,2,0,0))
parts.append(part11[part1[1]]); boolean_list.append(1)
#savestl('t1',part11,part1[1])
### part2
profile1 = Load('Text coordinates/t_prf_1','profile')
profile2 = Load('Text coordinates/t_prf_2','profile')
path = Load('Text coordinates/tip2', 'path')
part2 = sweep(np.array(profile1),np.array(profile2),path,0)
part22 = np.array(Transform(part2[0],0,180,0,-56,0,45))
parts.append(part22[part2[1]]); boolean_list.append(1)
#savestl('t2',part22,part2[1])
### part3
path3 = Load('Text coordinates/t_path3','path')
profile1 = [[2*math.cos(t/50),3.5*math.sin(t/50)] for t in range(int(2*math.pi*50))]
profile1.append(profile1[0])
part3 = sweep(np.array(profile1),np.array(profile1),path3,0)
part33 = np.array(Transform(part3[0],-90,90,0,28,0,45))
parts.append(part33[part3[1]]); boolean_list.append(1)
#savestl('t3',part33,part3[1])

