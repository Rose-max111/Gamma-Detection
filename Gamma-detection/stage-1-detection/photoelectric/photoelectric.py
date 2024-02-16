import numpy as np
import random
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde
import sys
import argparse


def init_cs(Z):
    file_path = f"./data/photoelectric/pe-cs-{Z}.dat"
    with open(file_path, "r") as file:
        lines = file.readlines()
    new_lines = [[float(num) for num in line.strip().split(" ")] for line in lines]
    Lower_bound, Upper_bound = new_lines[0][0], new_lines[0][1]

    ret_array = np.array(new_lines[2:])
    return Lower_bound, Upper_bound, ret_array


def init_le_cs(Z):
    file_path = f"./data/photoelectric/pe-le-cs-{Z}.dat"
    with open(file_path, "r") as file:
        lines = file.readlines()
    if len(lines) == 0:  # 避免 H 原子空文件的影响
        return 0, 0, 0
    new_lines = [[float(num) for num in line.strip().split(" ")] for line in lines]
    Lower_bound, Upper_bound = new_lines[0][0], new_lines[0][1]

    ret_array = np.array(new_lines[2:])
    return Lower_bound, Upper_bound, ret_array


def init_high(Z):
    file_path = f"./data/photoelectric/pe-high-{Z}.dat"
    with open(file_path, "r") as file:
        lines = file.readlines()
    new_lines = [[float(num) for num in line.strip().split(" ")] for line in lines]

    total_shell, use_shell, EA = new_lines[0][0], new_lines[0][1], new_lines[0][2]

    ret_array = np.array(new_lines[1:])
    return total_shell, use_shell, EA, ret_array


def init_low(Z):
    file_path = f"./data/photoelectric/pe-low-{Z}.dat"
    with open(file_path, "r") as file:
        lines = file.readlines()
    new_lines = [[float(num) for num in line.strip().split(" ")] for line in lines]

    total_shell, use_shell, EB = new_lines[0][0], new_lines[0][1], new_lines[0][2]

    ret_array = np.array(new_lines[1:])
    return total_shell, use_shell, EB, ret_array


class photoelectric:
    def __init__(self, Z):
        """
        储存所有的输入数据, 按缩写存储在相应变量中
        """
        self.cs_upper_bound, self.cs_lower_bound, self.cs_data = init_cs(Z)
        self.le_cs_upper_bound, self.le_cs_lower_bound, self.le_cs_data = init_le_cs(Z)
        (
            self.high_total_shell,
            self.high_used_shell,
            self.high_EA,
            self.high_data,
        ) = init_high(Z)
        (
            self.low_total_shell,
            self.low_used_shell,
            self.low_EB,
            self.low_data,
        ) = init_low(Z)

    def interpolation(self, E_gamma, data):
        """
        传入光子能量E_gamma, 插值数据集data

        返回该光子能量对应的散射截面大小(插值得到)
        """
        index = np.searchsorted(data[:, 0], E_gamma)
        if E_gamma == data[index, 0]:
            return data[index, 1] / (E_gamma**3)
        # data[:,0]有序，找到光子能量在其中的前驱后继进行插值
        xpre, ypre = data[index - 1, 0], data[index - 1, 1]
        xsuf, ysuf = data[index, 0], data[index, 1]
        return (ysuf - (ysuf - ypre) / (xsuf - xpre) * (xsuf - E_gamma)) / (
            E_gamma**3
        )

    def get_single_atom_photoelectric_cross_section(self, E_gamma):
        """
        传入光子能量E_gamma

        返回该光子能量对应的散射截面大小(根据不同光子能量分为若干区间进行了处理)
        """
        if E_gamma > self.high_EA:
            E_array = np.array([1 / (E_gamma ** (i + 1)) for i in range(6)])
            return np.sum(E_array * self.high_data[-1, 1::])
        elif E_gamma > self.low_EB:
            E_array = np.array([1 / (E_gamma ** (i + 1)) for i in range(6)])
            return np.sum(E_array * self.low_data[-1, 1::])
        elif E_gamma > self.high_data[0, 0]:
            return self.interpolation(E_gamma, self.cs_data)
        else:
            return self.interpolation(E_gamma, self.le_cs_data)

    def calculate_g(self, A, beta, gamma, nve):
        """
        传入A, beta, gamma, nve

        返回g(nve)的计算值
        """
        return (2 - nve) * (
            1 / (A + nve) + 1 / 2 * beta * gamma * (gamma - 1) * (gamma - 2)
        )

    def sampleing(self, E_gamma):
        """
        传入光子能量E_gamma

        利用特定抽样方法, 在满足该光电子角分布的情况下, 随机抽样出光电子的出射角theta

        返回该出射角的cos(theta) 以及 电子能量
        """
        E_elec = E_gamma - self.high_data[0, 0]
        m_e_c2 = 0.5109989461

        gamma = 1 + E_elec / m_e_c2
        beta = np.sqrt(E_elec * (E_elec + 2 * m_e_c2)) / (E_elec + m_e_c2)
        kappa = E_gamma / m_e_c2
        A = 1 / beta - 1
        g_0 = self.calculate_g(A, beta, gamma, 0)

        while True:
            kexi = random.random()
            nve = ((2 * A) / ((A + 2) ** 2 - 4 * kexi)) * (
                2 * kexi + (A + 2) * np.sqrt(kexi)
            )

            kexi_prime = random.random()
            g_nve = self.calculate_g(A, beta, gamma, nve)
            if kexi_prime * g_0 < g_nve:
                return 1 - nve, E_elec
        None


def draw_photoelectric_theta(Z, E_gamma):
    """
    传入原子核电荷大小Z(Z in 1/6/7/8/16)

    传入光子能量E_gamma

    作出该原子核, 光子能量在E_gamma时对应的光电子角分布曲线,并储存在pe_theta_{Z}.pdf中
    """
    assert Z in [1, 6, 7, 8, 16]
    atom = photoelectric(Z)
    output = np.array(
        [atom.sampleing(E_gamma)[0] for x in range(20000)]
    )  # 抽样20000次 以方便后面进行平滑处理
    output = np.arccos(output)
    output = output * 180 / np.pi

    # 以下使用了核密度估计 对离散的抽样光电子角分布数据进行平滑处理
    kde = gaussian_kde(output)
    x_values = np.linspace(0, 180, 2000)
    # 为了保证立体角，还需要再乘sin(theta)
    smoothed_density = kde(x_values) * np.sin(x_values / 180 * np.pi)
    # 归一化处理
    area_under_curve = np.trapz(smoothed_density, x=x_values)
    normalized_density = smoothed_density / area_under_curve

    atom_name = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S"}
    plt.plot(x_values, normalized_density, label=atom_name[Z])


def draw_all(E_gamma):
    for Z in [1, 6, 7, 8, 16]:
        draw_photoelectric_theta(Z, E_gamma)
    plt.xlabel("Theta(Degree)")
    plt.ylabel("dN/N")
    plt.title(f"Electron Angular Distribution At {E_gamma} MeV")
    plt.legend()
    plt.grid(True)
    plt.savefig(figure_name + ".pdf", format="pdf")


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

    draw_all(E_gamma)
