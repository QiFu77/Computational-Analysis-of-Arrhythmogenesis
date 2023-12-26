#coding=utf-8
import sys

from matplotlib import pyplot as plt
from matplotlib.pyplot import MultipleLocator  
import numpy as np
    
 


 
filename_con = "VentOneDResults_CON.dat"
data_con = np.loadtxt(filename_con).T[1:,:]
time = np.loadtxt(filename_con).T[0,:]
print(time[2502])
print(data_con.shape) # 2001,100

sys.exit(0)
'''
filename_corm = "VentOneDResults_CON3.dat"
data_corm = np.loadtxt(filename_corm).T[1:,:]
print(data_corm.shape) # 2001,100
'''

step_con_num = len(data_con[0,:]) # time
cell_con_num = len(data_con[:,0])

space = np.arange(0, cell_con_num, 1)       # Space range with space step
print(space.shape)


'''plt.figure(0)
# plt.contourf(time,space,zi,n,cmap=plt.cm.jet)
im = plt.contourf(time,space,data_con,500,cmap=plt.cm.jet,vmin = -81, vmax = 40)
plt.clim(-81,40)
# plt.contourf(time,space,data,cmap=plt.cm.jet) 
plt.colorbar(im)        
# plt.show()
'''

plt.figure(1,figsize=(5,4))
'''
ax1 = plt.subplot(111)
ax1.set_xlim(0,1000)
ax1.set_ylim(-0.05,0.25)
'''
plt.xlim(0,1000)
plt.ylim(-0.05,0.25)
#ax1.yaxis.set_major_locator(MultipleLocator(0.1))
plt.ylabel("Φ (mV)")
# dx = 0.15e-3 #m
# x0 = 20e-3 #m
# alpha = 11e-6 #m

dx = 0.15 #mm
x0 = 20+15 #mm
alpha = 11e-3 #mm

ECGx = []
ECG_cony = []
#ECG_cormy = []



for t in range(0,step_con_num,1):
	fai_con = 0
	for i in range(0,cell_con_num,1):

		if(i == 0):
			fai_con += -(data_con[i+1,t] - data_con[i,t])/((i*dx-x0)*(i*dx-x0))
		#	fai_corm += -(data_corm[i+1,t] - data_corm[i,t])/((i*dx-x0)*(i*dx-x0))

		elif(i == cell_con_num - 1):
			fai_con += -(data_con[i,t] - data_con[i-1,t])/((i*dx-x0)*(i*dx-x0))
			#fai_corm += -(data_corm[i,t] - data_corm[i-1,t])/((i*dx-x0)*(i*dx-x0))
			'''
			fai_INa += -(data_INa[i,t] - data_INa[i-1,t])/((i*dx-x0)*(i*dx-x0))
			fai_INaL += -(data_INaL[i,t] - data_INaL[i-1,t])/((i*dx-x0)*(i*dx-x0))
			fai_ICaL += -(data_ICaL[i,t] - data_ICaL[i-1,t])/((i*dx-x0)*(i*dx-x0))
			fai_IKr += -(data_IKr[i,t] - data_IKr[i-1,t])/((i*dx-x0)*(i*dx-x0))
			fai_IK1 += -(data_IK1[i,t] - data_IK1[i-1,t])/((i*dx-x0)*(i*dx-x0))
			'''
		else:
			fai_con += -(data_con[i+1,t] - data_con[i-1,t])/(2*(i*dx-x0)*(i*dx-x0))
			#fai_corm += -(data_corm[i+1,t] - data_corm[i-1,t])/(2*(i*dx-x0)*(i*dx-x0))
			'''
			fai_INa += -(data_INa[i+1,t] - data_INa[i-1,t])/(2*(i*dx-x0)*(i*dx-x0))
			fai_INaL += -(data_INaL[i+1,t] - data_INaL[i-1,t])/(2*(i*dx-x0)*(i*dx-x0))
			fai_ICaL += -(data_ICaL[i+1,t] - data_ICaL[i-1,t])/(2*(i*dx-x0)*(i*dx-x0))
			fai_IKr += -(data_IKr[i+1,t] - data_IKr[i-1,t])/(2*(i*dx-x0)*(i*dx-x0))
			fai_IK1 += -(data_IK1[i+1,t] - data_IK1[i-1,t])/(2*(i*dx-x0)*(i*dx-x0))
			'''
		# fai = fai*alpha*alpha/4
	ECGx.append(time[t]-7200)
	ECG_cony.append(fai_con)
	#ECG_cormy.append(fai_corm)


plt.xlim(0,1000)
plt.plot(ECGx,ECG_cony,color='black',linestyle='-',label='CON')
#ax1.plot(ECGx,ECG_cormy,color='red',linestyle='-',label='CORM-2',lw=1.3)

'''ax1.plot(ECGx,ECG_INay,color='MediumVioletRed',linestyle='-',label='$I\mathrm{_{Na}}$',LW = 1.3)
ax1.plot(ECGx,ECG_INaLy,color='SteelBlue',linestyle='-',label='$I\mathrm{_{NaL}}$',LW = 1.3)
ax1.plot(ECGx,ECG_ICaLy,color='SeaGreen',linestyle='-',label='$I\mathrm{_{CaL}}$',LW = 1.3)
ax1.plot(ECGx,ECG_IK1y,color='DarkOrange',linestyle='-',label='$I\mathrm{_{K1}}$',LW = 1.3)
ax1.plot(ECGx,ECG_IKry,color='Blue',linestyle='-',label='$I\mathrm{_{Kr}}$',LW = 1.3)
'''


#ax1.set_xticks([])
#plt.xticks([])#关闭坐标轴
# set range
# ax1.set_yticks(range(-0.2,0.31,0.1))
plt.ylim(-0.05,0.3)
#ax1.set_xticks(range(0,801,200))
# change label fontsize
#labels = ax1.get_xticklabels()+ax1.get_yticklabels()
#[label.set_fontsize(13) for label in labels]


'''ax2 = plt.subplot(212)
ax2.set_xlim(0,800)
#ax2.plot(ECGx,data_con0)
ax2.plot(ECGx,data_con0,'-.',color='black',label='ENDO')
ax2.plot(ECGx,data_con40,'-',color='black',label='MCELL')
ax2.plot(ECGx,data_con99,'--',color='black',label='EPI')
ax2.spines['top'].set_color('none')
ax2.spines['right'].set_color('none')
# set range
ax2.set_xticks(range(0,801,200))
ax2.set_yticks(range(-90,41,40))
# change label fontsize
labels = ax2.get_xticklabels()+ax2.get_yticklabels()
[label.set_fontsize(13) for label in labels]
'''

plt.legend(fontsize = 15, loc = 'best',frameon=False)
#plt.tight_layout()
plt.savefig('SQT6_heter.png',dpi=600)
plt.show()