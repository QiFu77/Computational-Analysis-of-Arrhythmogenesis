import numpy as np
import matplotlib.pyplot as plt

def plot_voltage_curves(files_and_params):
    plt.figure(figsize=(10, 6))
    for file_name, line_color, line_style, label_name, label_size, label_font in files_and_params:
        data = np.loadtxt(file_name)  # 不需要转置了
        time = data[:, 0]  # 时间在第一列
        voltage = data[:, 1]  # 膜电压在第二列

        plt.plot(time, voltage, color=line_color, linestyle=line_style, linewidth= 2.5, label=label_name,)  #设置线条粗细

    plt.ylabel('Membrane Potential (mV)',fontsize=20, fontweight='bold',fontproperties='Times New Roman')        #设置y轴标签的
    plt.yticks(np.arange(-100, 80, 40))
    plt.yticks(fontproperties='Times New Roman', size=20, weight='bold')  # 设置大小及加粗,y刻度轴
    plt.tick_params(axis='y', which='both', width=2)  # 设置y轴刻度线宽度
    plt.title('')
    plt.grid(False)  # 不显示网格线
    plt.gca().spines['top'].set_visible(False)  # 隐藏上侧边框
    plt.gca().spines['right'].set_visible(False)  # 隐藏右侧边框
    plt.gca().spines['bottom'].set_visible(False)  # 隐藏下侧边框
    plt.gca().spines['left'].set_linewidth(2)  # 设置y轴粗细
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)  # 隐藏x轴及其刻度值
    # 设置图例字体属性
    legend_properties = {'family': 'Times New Roman', 'weight': 'bold', 'size': 20}  #size 20 or 16
    plt.legend(loc='best',frameon=False, prop=legend_properties)  # 显示图例
    # plt.legend(frameon=False,fontsize=12)  # 显示图例
    plt.xlim(0, 800)
    #plt.savefig('TP06-wt+hete-Ab30&60-MCELL2',dpi=800)  #保存图片，名称根据要求修改
    #plt.savefig('TPORd-line2-with-Ab60&300-MCELL', dpi=800)  # 保存图片，名称根据要求修改
    plt.show()

files_and_params = [
    # ("TPORd/TPORd-HOMO-膜电位/VentriSingleCellResults_TPORd_HOMO_EPI.dat", 'red', '-', 'Homozygous', 20, 'Times New Roman'),
    ("TPORd/TPORd-WT-单细胞数据膜电压和erp/VentriSingleCellResults_TPORd_WT_MCELL.dat", 'grey', '--', 'Wild-type', 20, 'Times New Roman'),
    # ("TPORd/TPORd-WT-单细胞数据膜电压和erp/VentriSingleCellResults_TPORd_WT_MCELL.dat", 'black', '-', 'Wild-type', 20, 'Times New Roman'),
    ("TPORd/TPORd-dominate-膜电位/VentriSingleCellResults_TPORd_dominate_MCELL.dat", 'orange', '-', 'Heterozygous', 20, 'Times New Roman'),
    ("TPORd/TPORd-antibody60-膜电位+APD/VentriSingleCellResults_TPORd_dominate_MCELL_antibody60.dat", 'green', '-', 'KCNQ1 Ab 60μg/mg', 20, 'Times New Roman'),
    ("TPORd/TPORd-antibody300-膜电位-APD-ERP/VentriSingleCellResults_TPORd_dominate_MCELL_antibody300.dat", 'blue', '-', 'KCNQ1 Ab 300μg/mg', 20, 'Times New Roman'),
    #("VentOneDResultss_TPORd_different-BCL_DOMINATE_700.dat", 'blue', '-', 'DOMINATE-TPORd-700', 14, 'Times New Roman'),

    #-----------------------------------------TP06
    # ("TP06/TP06-全类型膜电位-APD/VentriSingleCellResults_WT_MCELL.dat", 'grey', '--', 'Wild-type', 15, 'Times New Roman'),
    #("TP06/TP06-全类型膜电位-APD/VentriSingleCellResults_WT_MCELL.dat", 'black', '-', 'Wild-type', 15, 'Times New Roman'),
    # ("TP06/TP06-全类型膜电位-APD/VentriSingleCellResults_dominate_MCELL.dat", 'orange', '-', 'Heterozygous', 15, 'Times New Roman'),
    # ("TP06/TP06-全类型膜电位-APD/VentriSingleCellResults_dominate_MCELL_antibody30.dat", 'green', '-', 'KCNQ1 Ab 30μg/mg', 15, 'Times New Roman'),
    # ("TP06/TP06-全类型膜电位-APD/VentriSingleCellResults_dominate_MCELL_antibody60.dat", 'blue', '-', 'KCNQ1 Ab 60μg/mg', 15, 'Times New Roman'),
    # ("TP06/TP06-全类型膜电位-APD/VentriSingleCellResults_HOMO_EPI.dat", 'red', '-', 'Homozygous', 15, 'Times New Roman'),
    # ("TP06/TP06-全类型膜电位-APD/0%移动-10膜电压EPI.dat", 'green', '-', 'KCNQ1 Ab 30μg/mg', 15, 'Times New Roman'),
    # ("TP06/TP06-全类型膜电位-APD/0%移动-28.1膜电压EPI.dat", 'blue', '-', 'KCNQ1 Ab 60μg/mg', 15, 'Times New Roman'),
    # 添加更多文件和参数
]

plot_voltage_curves(files_and_params)
