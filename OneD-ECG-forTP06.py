from matplotlib import pyplot as plt
from matplotlib.pyplot import MultipleLocator
import numpy as np

def plot_ECG(filename, BCL, label, color, linestyle, cycle_count):
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

        ECGx.append(time[t] - cycle_count * BCL)  # Use the cycle_count parameter
        ECGy.append(fai)

    return ECGx, ECGy

# File names and parameters
# Each tuple now includes an additional parameter for the cycle count
files_and_params = [
    # ("TPORd/TPORd-WT-differentBCL/VentOneDResultss_TPORd_different-BCL_WT_1000.dat", 1000, 'Wild-type', 'grey', '--', 59),
    # ("TPORd/TPORd-dominate-differentBCL/VentOneDResultss_TPORd_different-BCL_DOMINATE_1000.dat", 1000, 'Heterozygous', 'orange', '-', 59),
    # ("TPORd/TPORd-antibody60-ECG/VentOneDResultss_TPORd_different-BCL_DOMINATE_antibody60_1000.dat", 1000, 'KCNQ1 Ab 60μg/mg', 'green', '-', 59),
    # ("TP06/TP06-ECG-三种细胞/VentOneDResults_CON.dat", 1000, 'Wild-type', 'grey', '--', 10),
    ("TP06/TP06-ECG-三种细胞/VentOneDResults_CON.dat", 1000, 'Wild-type', 'grey', '--', 10),
    ("TP06/TP06-ECG-三种细胞/VentOneDResults_DOMINATE.dat", 1000, 'Heterozygous', 'orange', '-', 1),
    ("TP06/TP06-ECG-三种细胞/VentOneDResults_DOMINATE_antibody30.dat", 1000, 'KCNQ1 Ab 30μg/mg', 'green', '-', 1),
    ("TP06/TP06-ECG-三种细胞/VentOneDResults_DOMINATE_antibody60.dat", 1000, 'KCNQ1 Ab 60μg/mg', 'blue', '-', 1),
    # ("TP06/TP06-ECG-三种细胞/VentOneDResults_CON0%.dat", 1000, 'Homozygous', 'red', '-', 10),
    # ("TP06/TP06-ECG-三种细胞/0%移动-10mv.dat", 1000, 'KCNQ1 Ab 30μg/mg', 'green', '-', 10),
    # ("TP06/TP06-ECG-三种细胞/0%移动-28.1mv.dat", 1000, 'KCNQ1 Ab 60μg/mg', 'blue', '-', 10),

]

plt.figure(1, figsize=(10, 6))  # Adjust figure size to make it wider

ax = plt.gca()
ax.yaxis.set_major_locator(MultipleLocator(0.1))
plt.ylim((-0.05, 0.25))
plt.ylabel("Φ (mV)", fontsize=15, fontweight='bold', fontproperties='Times New Roman')  # 设置y轴标签

# Plot each ECG
for filename, BCL, label, color, linestyle, cycle_count in files_and_params:
    ECGx, ECGy = plot_ECG(filename, BCL, label, color, linestyle, cycle_count)
    plt.plot(ECGx, ECGy, color=color, label=label, linewidth=1.5, linestyle=linestyle)  # 增加了线条线条样式的功能

plt.xticks([])  # Disable x-axis ticks
plt.yticks(fontsize=15)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.ylim(-0.05, 0.3)
plt.xlim(0, 800)  # x轴的范围根据你的最大周期来确定

legend_properties = {'family': 'Times New Roman', 'weight': 'bold', 'size': 12}
plt.legend(loc='best', frameon=False, prop=legend_properties)  # 显示图例
plt.yticks(fontproperties='Times New Roman', size=15, weight='bold')  # 设置大小及加粗,y刻度轴
plt.tight_layout()
plt.savefig('TP06-ECG-WT+hete-Ab30&60', dpi=800)  # 保存图片，名称根据要求修改
plt.show()
