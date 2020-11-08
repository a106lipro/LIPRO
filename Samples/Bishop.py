# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 11:00:41 2020

@author: Vahid
"""

from LIPRO import *

profile1 = Load('text coordinates/Bishop','profile')
path = polygon(4,50,90)
path = Transform(np.insert(path,2,0,axis=1),0,90,0,-4,0,0)
part1 = sweep(profile1,profile1,path)
