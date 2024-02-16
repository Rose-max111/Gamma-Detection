from scipy import interpolate
from scipy.constants import c, electron_mass
import random
import math
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde
import os
import argparse

energy = electron_mass * c**2 / 1.6e-13  # 能量单位：MeV


def compton_sigma_e(E0, Z):
    """
    计算康普顿散射截面：

    输入: 入射光子能量（单个浮点数）、质子数（必须是整数）
    """
    with open("./data/compton/ce-cs-{}.dat".format(Z), "r") as f:
        first_line = f.readline().strip().split()
        Eu = float(first_line[1])
        El = float(first_line[0])
        data = f.readlines()[2:]
        x, y = zip(*[line.strip().split() for line in data])
        x = list(map(float, x))
        y = list(map(float, y))
    # 创建插值函数
    f = interpolate.interp1d(x, y)
    if E0 < 0.0001:
        return 0
    elif E0 < El:
        return E0 / El * compton_sigma_e(El)
    elif E0 > Eu:
        return Eu / E0 * compton_sigma_e(Eu)
    else:
        return f(E0) / E0

def deposition_point(E0, theta):
    """
    计算反冲电子出射角

    输入: 入射光子能量E0(单个浮点数)、 散射光子散射角度theta

    输出: theta_e: 电子出射角
    """
    theta_e = math.atan(1 / ((1 + E0 / energy) * math.tan(theta / 2)))
    return theta_e

def compton_scattered_photon_sampling(E0):
    """
    进行散射光子能量采样

    输入: 入射光子能量（单个浮点数)

    输出: 小epsilon(E1/E0) ; theta: 光子出射角; theta_e: 电子出射角
    """
    epsilon_0 = energy / (energy + 2 * E0)
    while True:
        r = random.random()
        alpha_1 = math.log(1 / epsilon_0)
        alpha_2 = 1 / 2 * (1 - epsilon_0**2)
        if r < alpha_1 / (alpha_1 + alpha_2):
            choice = 0
        else:
            choice = 1
        r_prime = random.random()
        if choice:
            epsilon = math.sqrt(epsilon_0**2 + (1 - epsilon_0**2) * r_prime)
        else:
            epsilon = epsilon_0**r_prime
        t = energy * (1 - epsilon) / (E0 * epsilon)
        sinsquare = t * (2 - t)
        cos = 1 - t
        theta = math.acos(cos)
        r_prime_prime = random.random()
        g = 1 - epsilon * sinsquare / (1 + epsilon**2)
        if g >= r_prime_prime:
            break
    theta_e = deposition_point(E0, theta)
    return epsilon, theta, theta_e


def compton_scattered_photon_ploting(E0):
    """
    绘制能量分布和角分布
    """
    energy_1_sample = []
    theta_sample = []
    theta_e_sample = []
    for i in range(1, 5000):
        energy_1_sample.append(compton_scattered_photon_sampling(E0)[0])
        theta_sample.append(compton_scattered_photon_sampling(E0)[1])
        theta_e_sample.append(compton_scattered_photon_sampling(E0)[2])
    energy_1 = np.array(energy_1_sample) * E0
    theta = np.array(theta_sample)
    theta_e = np.array(theta_e_sample)

    # 创建一个包含三个子图的画布和坐标系
    fig, (ax1, ax2, ax3) = plt.subplots(3)

    energy_1_smooth = np.linspace(min(energy_1), max(energy_1), 1000)
    kde_e1 = gaussian_kde(energy_1)
    ax1.plot(energy_1_smooth, kde_e1(energy_1_smooth))
    ax1.set_title("Energy 1")
    ax1.set_xlabel("Energy")
    ax1.set_ylabel("Frequency")

    # 绘制第二幅图
    theta_smooth = np.linspace(0, np.pi, 1000)
    kde_theta = gaussian_kde(theta)
    ax2.plot(theta_smooth, kde_theta(theta_smooth))
    ax2.set_title("Theta")
    ax2.set_xlabel("Theta")
    ax2.set_ylabel("Frequency")

    # 绘制第三幅图
    theta_e_smooth = np.linspace(0, np.pi / 2, 1000)
    kde_theta_e = gaussian_kde(theta_e)
    ax3.plot(theta_e_smooth, kde_theta_e(theta_e_smooth))
    ax3.set_title("Theta_e")
    ax3.set_xlabel("Theta_e")
    ax3.set_ylabel("Frequency")

    # 设置子图之间的间距
    plt.tight_layout()

    # 保存图像
    plt.savefig(f"{figure_name}.pdf".format(E0))


if __name__ == "__main__":
    # 创建 ArgumentParser 对象
    parser = argparse.ArgumentParser()

    # 添加命令行参数
    parser.add_argument("-e", type=float)
    parser.add_argument("-o", type=str)

    # 解析命令行参数
    args = parser.parse_args()

    # 获取参数值
    E_gamma = args.e
    figure_name = args.o
    compton_scattered_photon_ploting(E_gamma)
