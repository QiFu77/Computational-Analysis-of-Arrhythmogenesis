#
# General Code Structure (GCS) for one-dimensional plot
# 
# A. Code sructure is re-organised using standard C++ to fit into my project.
# B. The friendly style makes the project easier understood.
# C. The version is more extendable whatever to further 1D,2D or another single cell model. 
# 
# Under Intellectual Property Protection.
# 
# 
# Author      : Shugang Zhang <zhangshugang@hotmail.com>
# Last update : 06-10-2018
# 

#coding=utf-8    
from matplotlib import pyplot as plt    
import numpy as np
import math as math

 
# ------------put your filename here for visulization------------

filename="TPORd/TPORd-dominate-VW/40/Transmural1D_s2@425.00.dat"
#filename="TP06/TP06-ECG-三种细胞/VentOneDResults_DOMINATE.dat"
	
# -------------read in data------------
time = np.loadtxt(filename)[:,0].T   # the first colum should be time in my code structure
# time = np.arange(0, step_num, 1)   # if above time colum is not available, uncomment this line to auto-generate x-axis
data = np.loadtxt(filename)[:,1:].T # the rest columns are membrane potential. T is transverse matrix for better visulization
single = data[99,:] # this is for  singlecell observation, first dimension is the cell index after above transverse


# -------------time information------------
step_num = len(time[:])
print("Total time = %.3f ms" %(time[len(time)-1]))


# -------------space information------------
cell_num = len(data[:,0])
space = np.arange(0, cell_num, 1)  # Space range with space step
print("Cell num =",space.shape[0])

# -------------1D visualization------------
plt.figure(0)
im = plt.contourf(time,space,data,50,cmap=plt.cm.jet)
plt.clim(-80,40)           #----------------
# plt.contourf(time,space,data,cmap=plt.cm.jet)
cbar=plt.colorbar(im,ticks=[-90,40])
#cbar=plt.colorbar(im,ticks=[-70, 40]) # plot the color bar
cbar.ax.set_yticklabels([])
cbar.outline.set_visible(False)
plt.xticks([])
plt.yticks([])
#plt.xlim(1000, 2000)            #--------------x轴限制
# 获取当前的Axes对象
ax = plt.gca()
# 去掉黑色边框
for spine in ax.spines.values():
    spine.set_visible(False)
# -------------Single cell AP observation------------
'''
plt.figure(1)
# for i in range(0,80):
# 	plt.plot(time,data[5*i,:]+10*i)
for i in range(0,50):
	plt.plot(time,data[i,:])
'''
plt.savefig("TPORd-dominate3-易感窗5",dpi=800)
#plt.savefig("TP06-hete-ECG的彩虹图",dpi=800)
plt.show()