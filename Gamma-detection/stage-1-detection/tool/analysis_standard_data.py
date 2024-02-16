from matplotlib import pyplot as plt    
from scipy.stats import gaussian_kde
import h5py
import numpy as np
from matplotlib.ticker import MaxNLocator
import plotly.graph_objects as go
import plotly.io as pio
import math
# from sklearn.preprocessing import MinMaxScaler
output_data = []
with h5py.File("deposit_train.h5", "r") as opt:
    dataset_names = opt["track_info"]
    points = []
    for collision in dataset_names:
        points.append([collision[2],collision[3],collision[4]])
    points_array = np.array(points) 
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
    pio.write_html(fig, 'deposit_standard.html')