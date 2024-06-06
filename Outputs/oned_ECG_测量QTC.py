from matplotlib import pyplot as plt
from matplotlib.pyplot import MultipleLocator
import numpy as np


def plot_ECG(filename, BCL, label, color):
    data = np.loadtxt(filename).T[1:, :]
    time = np.loadtxt(filename).T[0, :]
    step_con_num = len(data[0, :])  # time
    cell_con_num = len(data[:, 0])
    print(data.shape)
    ECGx = []
    ECGy = []

    dx = 0.15  # mm
    x0 = 20 + 15  # mm
    alpha = 11e-3  # mm

    for t in range(0, step_con_num, 1):
        fai = 0
        for i in range(0, cell_con_num, 1):
            if (i == 0):
                fai += -(data[i + 1, t] - data[i, t]) / ((i * dx - x0) * (i * dx - x0))
            elif (i == cell_con_num - 1):
                fai += -(data[i, t] - data[i - 1, t]) / ((i * dx - x0) * (i * dx - x0))
            else:
                fai += -(data[i + 1, t] - data[i - 1, t]) / (2 * (i * dx - x0) * (i * dx - x0))

        ECGx.append(time[t] - 59 * BCL)                 #**************59是动态变化的，因为我们是运行了60个周期，输出的是最后一个周期
        ECGy.append(fai)                                #**************即59实际上应该是numS1-1,numS1是运行的周期数

    # 寻找满足条件的斜率最小的点
    threshold = 200                #******************要找的是t波斜率最小的点处的切线，因此时间取得200ms以后
    min_slope = float('inf')
    min_slope_index = -1

    for i in range(len(ECGy) - 1):
        if ECGx[i] > threshold:
            dx = ECGx[i + 1] - ECGx[i]
            dy = ECGy[i + 1] - ECGy[i]
            slope = dy / dx

            if slope < min_slope:
                min_slope = slope
                min_slope_index = i

    x_intersect = None
    if min_slope_index != -1:
        x1 = ECGx[min_slope_index]
        y1 = ECGy[min_slope_index]
        slope = (ECGy[min_slope_index + 1] - ECGy[min_slope_index]) / (
                    ECGx[min_slope_index + 1] - ECGx[min_slope_index])
        intercept = y1 - slope * x1
        x_intersect = -intercept / slope

    return ECGx, ECGy, x_intersect


# File names and parameters             #*******这里的输入格式就是文件名，对应的BCL周期/ms  线条标签，颜色。
files_and_params = [

    # ("TPORd/TPORd-WT-differentBCL/VentOneDResultss_TPORd_different-BCL_WT_1000.dat", 1000, 'Wild-type', 'black'),
    #("TPORd/TPORd-dominate-differentBCL/VentOneDResultss_TPORd_different-BCL_DOMINATE_1000.dat", 1000, 'Heterozygous', 'orange'),
    # ("TPORd/TPORd-antibody60-ECG/VentOneDResultss_TPORd_different-BCL_DOMINATE_antibody60_1000.dat", 1000, '60', 'green'),
    ("TPORd/TPORd-antibody300-differentBCl（似乎没用）/VentOneDResultss_TPORd_different-BCL_DOMINATE_antibody300_1000.dat", 1000, '300', 'blue'),
    # ("VentOneDResultss_TPORd_different-BCL_DOMINATE_600.dat", 600, 'DOMINATE-TPORd-600', 'red'),
    # ("VentOneDResultss_TPORd_different-BCL_DOMINATE_700.dat", 700, 'DOMINATE-TPORd-700', 'red'),
    # ("VentOneDResultss_TPORd_different-BCL_DOMINATE_800.dat", 800, 'DOMINATE-TPORd-800', 'red'),
    # ("VentOneDResultss_TPORd_different-BCL_DOMINATE_900.dat", 900, 'DOMINATE-TPORd-900', 'red'),
    # ("VentOneDResultss_TPORd_different-BCL_DOMINATE_1000.dat", 1000, 'DOMINATE-TPORd-1000', 'red'),
    # ("VentOneDResultss_TPORd_different-BCL_DOMINATE_1100.dat", 1100, 'DOMINATE-TPORd-1100', 'red'),
    # ("VentOneDResultss_TPORd_different-BCL_DOMINATE_1200.dat", 1200, 'DOMINATE-TPORd-1200', 'red'),
    # ("VentOneDResultss_TPORd_different-BCL_DOMINATE_1300.dat", 1300, 'DOMINATE-TPORd-1300', 'red'),
    # ("VentOneDResultss_TPORd_different-BCL_DOMINATE_1400.dat", 1400, 'DOMINATE-TPORd-1400', 'red'),
    # ("VentOneDResultss_TPORd_different-BCL_DOMINATE_1500.dat", 1500, 'DOMINATE-TPORd-1500', 'red'),
]

plt.figure(1, figsize=(10, 6))  # Adjust figure size to make it wider

ax = plt.gca()
ax.yaxis.set_major_locator(MultipleLocator(0.1))
plt.ylim((-0.05, 0.25))
plt.ylabel("Φ (mV)", fontsize=20,fontweight='bold',fontproperties='Times New Roman')#设置y轴标签

# Plot each ECG
for filename, BCL, label, color in files_and_params:
    ECGx, ECGy, x_intersect = plot_ECG(filename, BCL, label, color)
    plt.plot(ECGx, ECGy, color=color, label=label, linewidth=2.5)

    if x_intersect is not None:
        plt.axvline(x=x_intersect, color=color, linestyle='--')
        print(f"x value at y=0 for {label} after the threshold:", x_intersect)
        QTc_values= x_intersect/np.sqrt(60/(60000/BCL))                   #*********************x_intersect就是测量的QT值
        print(f"QTc values for {label}:", QTc_values)           #*************用巴雷特公式进行校验后的QTc的值

plt.xticks([])  # Disable x-axis ticks
plt.yticks(fontsize=20)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.ylim(-0.05, 0.3)
plt.xlim(0, 1000)   #*************************************这个是输出最后一个周期所展示的最大长度，我们的BCL范围是600-1500因此1600合理
legend_properties = {'family': 'Times New Roman', 'weight': 'bold', 'size': 20}
plt.legend(loc='best',frameon=False, prop=legend_properties)  # 显示图例
plt.yticks(fontproperties='Times New Roman', size=20, weight='bold')  # 设置大小及加粗,y刻度轴
#plt.legend(fontsize=15, loc='best', frameon=False)
plt.tight_layout()
#plt.savefig('TPORd-ECG-max1000-hete',dpi=800)  #保存图片，名称根据要求修改
plt.show()
