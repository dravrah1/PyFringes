# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 12:50:18 2017

@author: harv
"""
import numpy as np
import math
from PyQt5 import QtCore, QtGui, uic, QtWidgets

def r(n,m,p):
    order=12
    temp=np.zeros(order+1)
    for s in range(0,(n-m)+1):
        pow = 2 * (n-s) - m
        #print s,pow
        R = (-1)**s * math.factorial(2*n - m - s) / (math.factorial(s)*math.factorial(n-s)*math.factorial(n-m-s))
        temp[pow] = temp[pow] + R
    return temp


def genZern(nterms):
    i=0
    order=12
    #temp = np.zeros((nterms,order))
    temp = []
    for n in range(nterms):
        for m in range(n,-1,-1):
            if m == 0 :                              
                temp.append([i,n,m,r(n,m,0),'none']) 
                i=i+1
            else:
                temp.append([i,n,m,r(n,m,0),'cos'])
                i=i+1
                temp.append([i,n,m,r(n,m,0),'sin'])
                i=i+1
            if i > nterms:
                break
        if i > nterms:
            break
    return temp
        
def calcZern(zerns, coeffs, nterms, p, theta) :
    summ=float(0)

    for i in range(nterms):
        if coeffs[i] != 0:
            term=zerns[i]
            n = term[1]
            m = term[2]
            arr = term[3]
            fctn = term[4]
            termsum=0
            for j in range(2*n-m+1):
                termsum = termsum + arr[j]*(p**j)
            if fctn=='sin':
                summ = summ + termsum * coeffs[i] * math.sin(m*theta)
            elif fctn=='cos':
                summ = summ + termsum * coeffs[i] * math.cos(m*theta)
            else:
                summ = summ + termsum * coeffs[i]
    return summ
    
def genZernArrays2(nterms,N,zerns):
    zernarrays = np.zeros((nterms,N,N))
    summ = float(0)
    for x in range(N):
        for y in range(N):
            xx=float(x)-N/2.0
            yy=float(y)-N/2.0
            p = ((xx/(N/2))**2 + (yy/(N/2))**2)**0.5
            theta = math.atan2(yy,xx)
            for i in range(nterms):
                term=zerns[i]
                n = term[1]
                m = term[2]
                arr = term[3]
                fctn = term[4]
                termsum = 0
                summ = 0
                for j in range(2*n-m+1):
                    termsum = termsum + arr[j]*(p**j)
                    if fctn=='sin':
                        summ = summ + termsum * 1 * math.sin(m*theta)
                    elif fctn=='cos':
                        summ = summ + termsum * 1 * math.cos(m*theta)
                    else:
                        summ = summ + termsum * 1           
                zernarrays[i,x,y]=summ
    return zernarrays
                
def genZernArrays(nterms,N,zerns,progressbar):
#    progressbar = QtWidgets.QProgressBar
    zernarrays = np.zeros((nterms,N,N))
    summ = float(0)
    coeffs0 = np.zeros(nterms)
    for x in range(N):
        for y in range(N):
            progressbar.setValue(int(round(float(x*N+y)/float(N*N)*100.0)))
            progressbar.update()
            xx=float(x)-N/2.0
            yy=float(y)-N/2.0
            p = ((xx/(N/2))**2 + (yy/(N/2))**2)**0.5
            theta = math.atan2(yy,xx)
            coeffs=coeffs0
            for i in range(nterms):
                coeffs=np.zeros(nterms)
                coeffs[i]=1
                summ = calcZern(zerns,coeffs,nterms,p,theta)
                zernarrays[i,x,y]=summ
    return zernarrays        
        
                


            