from matplotlib import pyplot as plt
from matplotlib.pyplot import MultipleLocator
import numpy as np

def plot_ECG(filename, BCL, label, color,linestyle):
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

        ECGx.append(time[t] - 59 * BCL)  #*******这个59是可以修改的，因我的测试文件的周期是60，输出的是最后一个周期
                                         #********因此在不同的代码中可以修改
        ECGy.append(fai)

    return ECGx, ECGy

# File names and parameters
#***************************************第一个变量是文件名，第二个是BCL的周期，第三个是标签名字，最后一个是颜色代码，第四个是线条样式
files_and_params = [

    # ("TPORd/TPORd-WT-differentBCL/VentOneDResultss_TPORd_different-BCL_WT_1000.dat", 1000, 'Wild-type', 'black'),
    #("TPORd/TPORd-dominate-differentBCL/VentOneDResultss_TPORd_different-BCL_DOMINATE_1000.dat", 1000, 'Heterozygous', 'orange'),

    ("TPORd/TPORd-WT-differentBCL/VentOneDResultss_TPORd_different-BCL_WT_1000.dat", 1000, 'Wild-type', 'grey','--'),
    # ("TPORd/TPORd-WT-differentBCL/VentOneDResultss_TPORd_different-BCL_WT_1000.dat", 1000, 'Wild-type', 'black', '-'),
    ("TPORd/TPORd-dominate-differentBCL/VentOneDResultss_TPORd_different-BCL_DOMINATE_1000.dat", 1000, 'Heterozygous', 'orange','-'),
   # ("TPORd/TPORd-antibody60-ECG/VentOneDResultss_TPORd_different-BCL_DOMINATE_antibody60_1000.dat", 1000, 'KCNQ1 Ab 60μg/mg', 'green','-'),
    ("TPORd/TPORd-antibody300-differentBCl（似乎没用）/VentOneDResultss_TPORd_different-BCL_DOMINATE_antibody300_1000.dat", 1000, 'KCNQ1 Ab 300μg/mg', 'blue','-'),
    # ("TPORd/TPORd-dominate-differentBCL/VentOneDResultss_TPORd_different-BCL_DOMINATE_600.dat", 600, 'DOMINATE-TPORd-600', 'red'),
    # ("TPORd/TPORd-dominate-differentBCL/VentOneDResultss_TPORd_different-BCL_DOMINATE_700.dat", 700, 'DOMINATE-TPORd-700', 'red'),
    # ("TPORd/TPORd-dominate-differentBCL/VentOneDResultss_TPORd_different-BCL_DOMINATE_800.dat", 800, 'DOMINATE-TPORd-800', 'red'),
    # ("TPORd/TPORd-dominate-differentBCL/VentOneDResultss_TPORd_different-BCL_DOMINATE_900.dat", 900, 'DOMINATE-TPORd-900', 'red'),
    # ("TPORd/TPORd-dominate-differentBCL/VentOneDResultss_TPORd_different-BCL_DOMINATE_1000.dat", 1000, 'DOMINATE-TPORd-1000', 'red'),
    # ("TPORd/TPORd-dominate-differentBCL/VentOneDResultss_TPORd_different-BCL_DOMINATE_1100.dat", 1100, 'DOMINATE-TPORd-1100', 'red'),
    # ("TPORd/TPORd-dominate-differentBCL/VentOneDResultss_TPORd_different-BCL_DOMINATE_1200.dat", 1200, 'DOMINATE-TPORd-1200', 'red'),
    # ("TPORd/TPORd-dominate-differentBCL/VentOneDResultss_TPORd_different-BCL_DOMINATE_1300.dat", 1300, 'DOMINATE-TPORd-1300', 'red'),
    # ("TPORd/TPORd-dominate-differentBCL/VentOneDResultss_TPORd_different-BCL_DOMINATE_1400.dat", 1400, 'DOMINATE-TPORd-1400', 'red'),
    # ("TPORd/TPORd-dominate-differentBCL/VentOneDResultss_TPORd_different-BCL_DOMINATE_1500.dat", 1500, 'DOMINATE-TPORd-1500', 'red'),
    # ("VentOneDResultss_ORd_different-BCL_DOMINATE_700.dat", 700, 'DOMINATE-700', 'k'),
    #("VentOneDResultss_TPO6_numS1=60_different-BCL_DOMINATE_600.dat", 600, 'DOMINATE-TP06-600', 'orange')

]

plt.figure(1, figsize=(10, 6))  # Adjust figure size to make it wider

ax = plt.gca()
ax.yaxis.set_major_locator(MultipleLocator(0.1))
plt.ylim((-0.05, 0.25))
plt.ylabel("Φ (mV)", fontsize=20,fontweight='bold',fontproperties='Times New Roman')#设置y轴标签

# Plot each ECG
for filename, BCL, label, color,linestyle in files_and_params:
    ECGx, ECGy = plot_ECG(filename, BCL, label, color,linestyle)
    plt.plot(ECGx, ECGy, color=color, label=label, linewidth=2.5,linestyle=linestyle)                     #增加了线条线条样式的功能

plt.xticks([])  # Disable x-axis ticks
plt.yticks(fontsize=20)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.ylim(-0.05, 0.3)
plt.xlim(0, 800)  #**********x轴的范围根据你的最大周期来确定

legend_properties = {'family': 'Times New Roman', 'weight': 'bold', 'size': 20}
plt.legend(loc='best',frameon=False, prop=legend_properties)  # 显示图例
plt.yticks(fontproperties='Times New Roman', size=20, weight='bold')  # 设置大小及加粗,y刻度轴
legend_properties = {'family': 'Times New Roman', 'weight': 'bold', 'size': 20}
#plt.legend(fontsize=15, loc='best', frameon=False)
plt.tight_layout()
plt.savefig('TPORd-ECG-WT-Ab300（虚线）',dpi=800)  #保存图片，名称根据要求修改
plt.show()