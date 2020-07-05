# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 16:53:42 2015

@author: Tony & Vic
"""

import matplotlib.pyplot as plt
from numpy import empty

ff = open('parameters1.txt', 'r')
gg = open('parameters2.txt', 'r')
hh = open('parameters3.txt', 'r')
#ff = open('telescope.txt', 'r')

########## retrieving r_200 and c from parameters.txt w/ endens = 30 ##########

c= [map(float, line.split('        ')) for line in ff ]

ff.close()

size = len(c)

x = empty([size, 1])
y = empty([size, 1])
for i in range(0,size):
    x[i] = c[i][0] # Save the first element of the ith galaxy
    y[i] = c[i][1] # Save the second element of the ith galaxy

########## retrieving r_200 and c from parameters2.txt w/ endens = 60 #########

d= [map(float, line.split('        ')) for line in gg ]

gg.close()

size1 = len(d)

x2 = empty([size1, 1])
y2 = empty([size1, 1])
for i in range(0,size1):
    x2[i] = d[i][0] # Save the first element of the ith galaxy
    y2[i] = d[i][1] # Save the second element of the ith galaxy   
    
######### retrieving r_200 and c from parameters3.txt w/ endens = 300 #########

e= [map(float, line.split('        ')) for line in hh ]

hh.close()

size1 = len(e)

x3 = empty([size1, 1])
y3 = empty([size1, 1])
for i in range(0,size1):
    x3[i] = e[i][0] # Save the first element of the ith galaxy
    y3[i] = e[i][1] # Save the second element of the ith galaxy   
        
    
########################## plotting them in 1 graph ###########################    
    
plt.plot(x,y,'ro')
plt.plot(x2,y2,'go')
plt.plot(x3,y3,'ko')
plt.legend(['Density of 30', 'Density of 60', 'Density of 300'])
plt.xlabel('r_200'); plt.ylabel('concentration')
