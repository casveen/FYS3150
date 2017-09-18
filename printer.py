# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 21:07:37 2017

@author: chris
"""
import matplotlib.pyplot as plt
import numpy as np

def fileToArray(filename):
     out=[]
     f=open(filename,'r');
     for line in f:
         out.append(float(line))
     return out


def plotDifferences(file1, file2):
    y1=fileToArray(file1)
    y2=fileToArray(file2)
    m=min(len(y1),len(y2))
    x=np.arange(0,1,1/m)
    
    plt.plot(x,y2,'b.')        
    plt.plot(x,y1,'r-')
    plt.legend(['exact','n=%d'%(len(y1))])
    plt.show()

def plot4Differences(file1, file2, file3, file4):
    y1=fileToArray(file1)
    y2=fileToArray(file2)
    y3=fileToArray(file3)
    y4=fileToArray(file4)
    
    plt.figure(figsize=(20,20))
    plt.plot(np.arange(0,1,1/len(y1)),y1,"r-")
    plt.plot(np.arange(0,1,1/len(y2)),y2,"b-")
    plt.plot(np.arange(0,1,1/len(y3)),y3,"y-")
    plt.plot(np.arange(0,1,1/len(y4)),y4,"g-")
    plt.legend(['n=10','n=100','n=1000','exact'])
    
    plt.show()