import numpy as np
import matplotlib.pyplot as plt

def plot_voltage_curves(files_and_params):
    plt.figure(figsize=(10, 6))
    for file_name, line_color, line_style, label_name, label_size, label_font in files_and_params:
        data = np.loadtxt(file_name)  # 不需要转置了
        time = data[:, 0]  # 时间在第一列
        voltage = data[:, 1]  # 膜电压在第二列

        plt.plot(time, voltage, color=line_color, linestyle=line_style, label=label_name)

    plt.ylabel('Membrane Voltage')
    plt.title('')
    plt.grid(False)  # 不显示网格线
    plt.gca().spines['top'].set_visible(False)  # 隐藏上侧边框
    plt.gca().spines['right'].set_visible(False)  # 隐藏右侧边框
    plt.gca().spines['bottom'].set_visible(False)  # 隐藏下侧边框
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)  # 隐藏x轴及其刻度值
    plt.legend(frameon=False)  # 显示图例
    plt.show()


files_and_params = [
    ("VentriSingleCellResults_TPORd_dominate_MCELL.dat", 'red', '-', 'ENDO', 12, 'Arial'),
   # ("VentOneDResultss_TPORd_different-BCL_DOMINATE_700.dat", 'blue', '-', 'DOMINATE-TPORd-700', 14, 'Times New Roman'),
    # 添加更多文件和参数
]

plot_voltage_curves(files_and_params)
