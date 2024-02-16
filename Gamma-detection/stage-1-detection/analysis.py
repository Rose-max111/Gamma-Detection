from matplotlib import pyplot as plt    
from scipy.stats import gaussian_kde
import h5py
import numpy as np
from matplotlib.ticker import MaxNLocator
import plotly.graph_objects as go
import plotly.io as pio
import math
import argparse
# from sklearn.preprocessing import MinMaxScaler
def draw_deposit(inputfile):
    output_data = []
    with h5py.File(inputfile, "r") as opt:
        group = opt["deposit_point"]
        dataset_names = group.keys()
        dataset_len = len(dataset_names)
        for t in range(dataset_len):
            tmp = opt["deposit_point"][f"{t}"]
            output_data.append(tmp)
        points = []
        dissymmetrypoints= []
        for photon in output_data:
            for collision in photon:
                    points.append([collision[1],collision[2],collision[3]])
                    dissymmetrypoints.append([collision[1],math.sqrt(collision[2]**2+collision[3]**2)])
        points_array = np.array(points) 
        dissymmetrypoints_array = np.array(dissymmetrypoints)
        # 从点数据数组中提取坐标
        x = points_array[:, 0]
        y = points_array[:, 1]
        z = points_array[:, 2]

        # 创建交互式3D散点图
        fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z, mode='markers', marker=dict(size=3, opacity=0.2))])

        # 设置图形布局
        fig.update_layout(scene=dict(
            xaxis_title='X 轴/mm',
            yaxis_title='Y 轴/mm',
            zaxis_title='Z 轴/mm',
            aspectmode='data'
        ))

        # 保存为HTML文件
        pio.write_html(fig, inputfile[:-3]+'.html')
        # 计算核密度估计
        kde = gaussian_kde(dissymmetrypoints_array.T)

        # 创建网格以绘制等高线图
        x_min, x_max = dissymmetrypoints_array[:, 0].min() - 0.1, dissymmetrypoints_array[:, 0].max() + 0.1
        y_min, y_max = dissymmetrypoints_array[:, 1].min() - 0.1, dissymmetrypoints_array[:, 1].max() + 0.1
        X, Y = np.mgrid[x_min:x_max:500j, y_min:y_max:500j]
        positions = np.vstack([X.ravel(), Y.ravel()])
        Z = np.reshape(kde(positions).T, X.shape)

        # 绘制等高线图
        plt.figure(figsize=(8, 6))
        levels = 100
        plt.contourf(X, Y, Z, levels, cmap='viridis')
        plt.colorbar(label='probability density')
        # plt.scatter(dissymmetrypoints_array[:, 0], dissymmetrypoints_array[:, 1], c='r', s=20, label='points')
        plt.xlabel('X/mm')
        plt.ylabel('distance/mm')
        plt.title('Probability density distribution on a two-dimensional plane')
        # plt.legend()
        plt.savefig(inputfile[:-3]+'_dissymmetry.png')
if __name__ == "__main__":
    # 创建 ArgumentParser 对象
    parser = argparse.ArgumentParser()

    # 添加命令行参数
    parser.add_argument("-i", type=str)
    # 解析命令行参数
    args = parser.parse_args()

    # 获取参数值
    inputfile = args.i  # 随机次数
    draw_deposit(inputfile)
