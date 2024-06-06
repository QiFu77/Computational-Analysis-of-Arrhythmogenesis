from matplotlib import pyplot as plt
from matplotlib.pyplot import MultipleLocator
import numpy as np

#是把一类个体的不同BCL的ECG画在一张图里，需要修改的将会是  1.filename_con,  2.每一个文件的BCL的值

#BCL1
filename_con =  "VentOneDResultss_TPORd_different-BCL_DOMINATE_600.dat"
data_con = np.loadtxt(filename_con).T[1:, :]
time = np.loadtxt(filename_con).T[0, :]
print(data_con.shape)  # 2001,100
step_con_num = len(data_con[0, :])  # time
cell_con_num = len(data_con[:, 0])
BCL1=600
ECGx = []
ECG_cony = []
#----------------------
#BCL2
filename_con2 =  "VentOneDResultss_ORd_different-BCL_DOMINATE_700.dat"
data_con2 = np.loadtxt(filename_con2).T[1:, :]
time2 = np.loadtxt(filename_con2).T[0, :]
print(data_con2.shape)  # 2001,100
step_con_num2 = len(data_con2[0, :])  # time
cell_con_num2 = len(data_con2[:, 0])
BCL2=700
ECGx2 = []
ECG_cony2 = []
#----------------------
#BCL3
filename_con3 =  "VentOneDResultss_ORd_different-BCL_DOMINATE_800.dat"
data_con3 = np.loadtxt(filename_con3).T[1:, :]
time3 = np.loadtxt(filename_con3).T[0, :]
print(data_con3.shape)  # 2001,100
step_con_num3 = len(data_con3[0, :])  # time
cell_con_num3 = len(data_con3[:, 0])
BCL3=800
ECGx3 = []
ECG_cony3 = []
#----------------------
#BCL4
filename_con4 =  "VentOneDResultss_ORd_different-BCL_DOMINATE_900.dat"
data_con4 = np.loadtxt(filename_con4).T[1:, :]
time4 = np.loadtxt(filename_con4).T[0, :]
print(data_con4.shape)  # 2001,100
step_con_num4 = len(data_con4[0, :])  # time
cell_con_num4 = len(data_con4[:, 0])
BCL4=900
ECGx4 = []
ECG_cony4 = []
#----------------------
#BCL5
filename_con5 =  "VentOneDResultss_ORd_different-BCL_DOMINATE_1000.dat"
data_con5 = np.loadtxt(filename_con5).T[1:, :]
time5 = np.loadtxt(filename_con5).T[0, :]
print(data_con5.shape)  # 2001,100
step_con_num5 = len(data_con5[0, :])  # time
cell_con_num5 = len(data_con5[:, 0])
BCL5=1000
ECGx5 = []
ECG_cony5 = []
#----------------------
#BCL6
filename_con6 =  "VentOneDResultss_ORd_different-BCL_DOMINATE_1100.dat"
data_con6 = np.loadtxt(filename_con6).T[1:, :]
time6 = np.loadtxt(filename_con6).T[0, :]
print(data_con6.shape)  # 2001,100
step_con_num6 = len(data_con6[0, :])  # time
cell_con_num6 = len(data_con6[:, 0])
BCL6=1100
ECGx6 = []
ECG_cony6 = []
#----------------------
#BCL7
filename_con7 =  "VentOneDResultss_TPORd_different-BCL_DOMINATE_1200.dat"
data_con7 = np.loadtxt(filename_con7).T[1:, :]
time7 = np.loadtxt(filename_con7).T[0, :]
print(data_con7.shape)  # 2001,100
step_con_num7 = len(data_con7[0, :])  # time
cell_con_num7 = len(data_con7[:, 0])
BCL7=1200
ECGx7 = []
ECG_cony7 = []
#----------------------


#BCL....
#----------------------

space = np.arange(0, cell_con_num, 1)  # Space range with space step
print(space.shape)

plt.figure(1, figsize=(6, 4))  # 调整图形尺寸，使其更宽

ax = plt.gca()
ax.yaxis.set_major_locator(MultipleLocator(0.1))
plt.ylim((-0.05, 0.25))
plt.ylabel("Φ (mV)", fontsize=15)

dx = 0.15  # mm
x0 = 20 + 15  # mm
alpha = 11e-3  # mm

#对每一个不同的BCL都重复的画图
#-----------------------------BCL1
for t in range(0, step_con_num, 1):
    fai_con = 0
    for i in range(0, cell_con_num, 1):
        if (i == 0):
            fai_con += -(data_con[i + 1, t] - data_con[i, t]) / ((i * dx - x0) * (i * dx - x0))
        elif (i == cell_con_num - 1):
            fai_con += -(data_con[i, t] - data_con[i - 1, t]) / ((i * dx - x0) * (i * dx - x0))
        else:
            fai_con += -(data_con[i + 1, t] - data_con[i - 1, t]) / (2 * (i * dx - x0) * (i * dx - x0))
    ECGx.append(time[t]-59*BCL1)
    ECG_cony.append(fai_con)
#-----------------------------BCL2
for t in range(0, step_con_num2, 1):
    fai_con2 = 0
    for i in range(0, cell_con_num2, 1):
        if (i == 0):
            fai_con2 += -(data_con2[i + 1, t] - data_con2[i, t]) / ((i * dx - x0) * (i * dx - x0))
        elif (i == cell_con_num2 - 1):
            fai_con2 += -(data_con2[i, t] - data_con2[i - 1, t]) / ((i * dx - x0) * (i * dx - x0))
        else:
            fai_con2 += -(data_con2[i + 1, t] - data_con2[i - 1, t]) / (2 * (i * dx - x0) * (i * dx - x0))
    ECGx2.append(time2[t]-59*BCL2)
    ECG_cony2.append(fai_con2)
#-----------------------------------------------
#-----------------------------BCL3
for t in range(0, step_con_num3, 1):
    fai_con3 = 0
    for i in range(0, cell_con_num3, 1):
        if (i == 0):
            fai_con3 += -(data_con3[i + 1, t] - data_con3[i, t]) / ((i * dx - x0) * (i * dx - x0))
        elif (i == cell_con_num3 - 1):
            fai_con3 += -(data_con3[i, t] - data_con3[i - 1, t]) / ((i * dx - x0) * (i * dx - x0))
        else:
            fai_con3 += -(data_con3[i + 1, t] - data_con3[i - 1, t]) / (2 * (i * dx - x0) * (i * dx - x0))
    ECGx3.append(time3[t]-59*BCL3)
    ECG_cony3.append(fai_con3)
#-----------------------------------------------
#-----------------------------BCL4
for t in range(0, step_con_num4, 1):
    fai_con4 = 0
    for i in range(0, cell_con_num4, 1):
        if (i == 0):
            fai_con4 += -(data_con4[i + 1, t] - data_con4[i, t]) / ((i * dx - x0) * (i * dx - x0))
        elif (i == cell_con_num4 - 1):
            fai_con4 += -(data_con4[i, t] - data_con4[i - 1, t]) / ((i * dx - x0) * (i * dx - x0))
        else:
            fai_con4 += -(data_con4[i + 1, t] - data_con4[i - 1, t]) / (2 * (i * dx - x0) * (i * dx - x0))
    ECGx4.append(time4[t]-59*BCL4)
    ECG_cony4.append(fai_con4)
#-----------------------------------------------
#-----------------------------BCL5
for t in range(0, step_con_num5, 1):
    fai_con5 = 0
    for i in range(0, cell_con_num5, 1):
        if (i == 0):
            fai_con5 += -(data_con5[i + 1, t] - data_con5[i, t]) / ((i * dx - x0) * (i * dx - x0))
        elif (i == cell_con_num5 - 1):
            fai_con5 += -(data_con5[i, t] - data_con5[i - 1, t]) / ((i * dx - x0) * (i * dx - x0))
        else:
            fai_con5 += -(data_con5[i + 1, t] - data_con5[i - 1, t]) / (2 * (i * dx - x0) * (i * dx - x0))
    ECGx5.append(time5[t]-59*BCL5)
    ECG_cony5.append(fai_con5)
#-----------------------------------------------
#-----------------------------BCL6
for t in range(0, step_con_num6, 1):
    fai_con6 = 0
    for i in range(0, cell_con_num6, 1):
        if (i == 0):
            fai_con6 += -(data_con6[i + 1, t] - data_con6[i, t]) / ((i * dx - x0) * (i * dx - x0))
        elif (i == cell_con_num6 - 1):
            fai_con6 += -(data_con6[i, t] - data_con6[i - 1, t]) / ((i * dx - x0) * (i * dx - x0))
        else:
            fai_con6 += -(data_con6[i + 1, t] - data_con6[i - 1, t]) / (2 * (i * dx - x0) * (i * dx - x0))
    ECGx6.append(time6[t]-59*BCL6)
    ECG_cony6.append(fai_con6)
#-----------------------------------------------
#-----------------------------BCL7
for t in range(0, step_con_num7, 1):
    fai_con7 = 0
    for i in range(0, cell_con_num7, 1):
        if (i == 0):
            fai_con7 += -(data_con7[i + 1, t] - data_con7[i, t]) / ((i * dx - x0) * (i * dx - x0))
        elif (i == cell_con_num7 - 1):
            fai_con7 += -(data_con7[i, t] - data_con7[i - 1, t]) / ((i * dx - x0) * (i * dx - x0))
        else:
            fai_con7 += -(data_con7[i + 1, t] - data_con7[i - 1, t]) / (2 * (i * dx - x0) * (i * dx - x0))
    ECGx7.append(time7[t]-59*BCL7)
    ECG_cony7.append(fai_con7)
#-----------------------------------------------

plt.plot(ECGx7, ECG_cony7, color='orange', label='DOMINATE-1200', linewidth=1.5)
plt.plot(ECGx6, ECG_cony6, color='silver', label='WT-1500', linewidth=1.5)
plt.plot(ECGx5, ECG_cony5, color='darkgray', label='1400', linewidth=1.5)
plt.plot(ECGx4, ECG_cony4, color='grey', label='1200', linewidth=1.5)
plt.plot(ECGx3, ECG_cony3, color='black', label='600', linewidth=1.5)
plt.plot(ECGx2, ECG_cony2, color='k', label='800', linewidth=1.5)
plt.plot(ECGx, ECG_cony, color='red', label='600', linewidth=1.5)

plt.xticks([])  # 关闭x轴坐标刻度
plt.yticks(fontsize=15)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.ylim(-0.05, 0.3)
plt.xlim(0,1600)#*******************右边界取BCL的最大值，本次测试不超过1500，因此取得是1600

plt.legend(fontsize=15, loc='best', frameon=False)
plt.tight_layout()
plt.show()
