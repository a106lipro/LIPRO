# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 20:38:49 2020

@author: Vahid
"""

from LIPRO import *

G1 = polygon(10,4)
S1 = [[0,0,0],[0,0,10*np.sqrt(2)]]
part1 = sweep(G1,G1,S1)

#savestl('cube',part1[0],part1[1])