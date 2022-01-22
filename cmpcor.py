#!/usr/bin/env python
#-*-coding:utf-8-*-
from obspy import read
import numpy as np
import matplotlib.pyplot as plt
import time
import copy
import sys

#if you have not obspy package,please take other method to gen data for oriData.


def Cor_Cal(npA,npB,begin):
     COV = np.dot(npA, npB)  
     return COV


def Cor_CalSeries(npA,npB):#npA领先的时间
   if len(npA)!=len(npB):
         print "Error len npA != len npB."
         return []
   lenX= len(npA)
   calUseTimeLen=np.zeros(lenX*2)
   y= np.zeros(lenX*3)
   y[lenX:lenX*2] = npA
   xmove= np.zeros(lenX*5)
   xmove[lenX*2:lenX*3] = npB
   for i in np.arange(lenX*2):
      x= xmove[2*lenX-i:5*lenX-i]#npB不断地前移，前移多少就是npA领先多少
      calUseTimeLen[i]=Cor_Cal(x,y,i)
   calUseTimeLen[0:len(calUseTimeLen)-1] = calUseTimeLen[1:len(calUseTimeLen)]
   return calUseTimeLen



#### step one prepare data #######
syspath = sys.path[0]
oriData=[]
if True:
  st = read(syspath+"/test1.ascii")
  for tr in st:
     oriData.append(tr.data) 
  st = read(syspath+"/test2.sac")
  for tr in st:
     oriData.append(tr.data)
if False:
  data1=np.random.random(2000)+np.sin(np.arange(2000)/20.0)*np.cos(np.arange(2000))
  data2=np.random.random(2000)+np.cos(np.arange(2000)/20.0)
  oriData.append(data1)
  oriData.append(data2)

if True:
  useData={}
  mean1=(oriData[0][0:1000]).mean()
  mean2=(oriData[1][0:1000]).mean()
  oriData[0]=(oriData[0][0:1000])-mean1
  oriData[1]=(oriData[1][0:1000])-mean2
  useData[0]= np.append(oriData[0],np.zeros(1000))
  useData[1]= np.append(oriData[1],np.zeros(1000))


####step two cal cor use fft#######
  tmp0 = np.fft.fft(useData[0])
  tmp1 = np.fft.fft(useData[1])
  tmp = tmp0*(tmp1.conjugate())#tmp0领先的时间
  out = np.fft.ifft(tmp)
  xt= np.arange(len(tmp))
  calUseFFT = np.append(out.real[1001:2000],out.real[0:1000])
####step three cal cor use time#########################

if True:
  calUseTime = Cor_CalSeries(oriData[0],oriData[1])
if True:
  if len(calUseFFT)>len(calUseTime):
     calUseFFT=calUseFFT[0:len(calUseTime)]
  else:
     calUseTime=calUseTime[0:len(calUseFFT)]
if True:
  delta1=(calUseFFT-calUseTime)
  print "time out Maxvalue:", max(abs(calUseTime))
  print "fft out Maxvalue:",max(abs(calUseFFT))
  print "two mode max bias between using fft and using Time:",delta1.max()

##############step four use python correlate #########

pythonOut=np.correlate(useData[0], useData[1],"same")#useData[0]领先的时间
pythonOut=pythonOut[1:]
delta2=(pythonOut-calUseTime)
print "two mode max bias between using np.correlate and using Time:",delta2.max()




####step five plot#######
fig = plt.figure(figsize=(25,14))
plt.subplots_adjust(wspace = 0.1,hspace =0.4,left=0.05,right=0.95)
ax11 = fig.add_subplot(4,2,1)
ax12 = fig.add_subplot(4,2,3)
ax13 = fig.add_subplot(4,2,5)
ax14 = fig.add_subplot(4,2,7)

ax21 = fig.add_subplot(4,2,2)
ax22 = fig.add_subplot(4,2,4)
ax23 = fig.add_subplot(4,2,6)
ax24 = fig.add_subplot(4,2,8)



ax11.plot(oriData[0])
ax11.set_title("original data A")


ax21.plot(oriData[1])
ax21.set_title("original data B")


ax12.plot(useData[0])
ax13.plot(useData[1])

lenX=len(oriData[0])
y= np.zeros(lenX*3)
y[lenX:lenX*2] =oriData[0]
ax22.plot(y)

lenX=len(oriData[1])
y= np.zeros(lenX*3)
y[lenX:lenX*2] =oriData[1]
ax23.plot(y)



xt= np.arange(len(calUseFFT))
xt = xt - len(calUseFFT)/2
ax14.plot(xt,calUseFFT)

xt= np.arange(len(calUseTime))
xt = xt - len(calUseFFT)/2
ax24.plot(xt,calUseTime)


ax12.set_title("basic data one")
ax13.set_title("basic data two")
ax14.set_title("basic data one-two cor use FFT")


ax22.set_title("basic data 1")
ax23.set_title("basic data 2")
ax24.set_title("basic data 1-2 use time cor")


ax12.grid()
ax13.grid()
ax14.grid()
ax22.grid()
ax23.grid()
ax24.grid()


lineX=np.arange(1000,2000)
lineY=np.zeros(1000)
ax12.plot(lineX,lineY,color="red",linewidth=2)
ax13.plot(lineX,lineY,color="red",linewidth=2)

lineX=np.arange(0,1000)
lineY=np.zeros(1000)
ax22.plot(lineX,lineY,color="red",linewidth=2)
ax23.plot(lineX,lineY,color="red",linewidth=2)


lineX=np.arange(2000,3000)
lineY=np.zeros(1000)
ax22.plot(lineX,lineY,color="red",linewidth=2)
ax23.plot(lineX,lineY,color="red",linewidth=2)
picsavename = syspath+"/cmpcor.png"
plt.savefig(picsavename)
print picsavename,"have be saved."
plt.show()
############################################





