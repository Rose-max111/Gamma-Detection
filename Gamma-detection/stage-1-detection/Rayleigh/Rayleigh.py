from collections import defaultdict
from scipy import interpolate
from scipy.constants import c, electron_mass, elementary_charge, epsilon_0, Planck
import random
import math
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde
import xml.etree.ElementTree as ET
import os
import argparse

re = (
    elementary_charge**2 / (4 * math.pi * epsilon_0 * electron_mass * c**2) * 100
)  # 电子经典半径，单位为cm
F = defaultdict(list)  # 质子数：F插值对象
avogadro_number = 6.022e23


def extract_density(filename):
    # 输入：文件名
    # 输出：一个字典。键值对为元素：数密度n
    tree = ET.parse(filename)
    root = tree.getroot()

    density_dict = {}
    n_value = []
    ref_name = []
    atom_value = []

    for fraction in root.iter("fraction"):
        n_value.append(float(fraction.attrib.get("n")))
        ref_name.append(fraction.attrib.get("ref"))

    D_value = float(root.find(".//material/D").get("value"))
    for atom in root.iter("atom"):
        atom_value.append(float(atom.attrib.get("value")))

    for i in range(5):
        density_dict[ref_name[i]] = (
            n_value[i] * D_value / atom_value[i] * avogadro_number
        )

    return density_dict


def Z_number(filename):
    # 输入：文件名
    # 输出：一个字典。键值对为元素：质子数Z
    tree = ET.parse(filename)
    root = tree.getroot()
    string = extract_density(filename).keys()
    Z_dict = {}

    for element in root.iter("element"):
        name = element.attrib.get("name")
        for prefix_length in range(1, len(name) + 1):
            prefix = name[:prefix_length]
            if prefix in string:
                break
        atom_value = int(element.attrib.get("Z"))
        Z_dict[prefix] = atom_value

    return Z_dict


def Rayleigh_sigma_e(E0, Z, F):
    """
    计算瑞利散射截面：

    输入: 入射光子能量（单个浮点数）、质子数（必须是整数）、插值函数F
    输出: 总瑞利散射截面
    """
    E0 = E0 * 1e6 * elementary_charge
    lamda = Planck * c / E0 * 1e10  # 波长单位: 埃
    F_x = F[Z]
    sum_integral = 0
    num_points = 100
    uniform_points = list(np.linspace(-1, 1, num_points))
    for cos_theta in uniform_points:
        sin_half_theta = math.sqrt((1 - cos_theta) / 2)
        x = sin_half_theta / lamda
        if x > 1e9:
            f = F_x(1e9) / 1e-18 * x**-2
        else:
            f = F_x(x)
        sum_integral += (
            0.02 * (1 + cos_theta**2) * f**2 * math.pi * re**2 * 1e24
        )  # 单位: barn
    return sum_integral


def get_Rayleigh_sigma(E_gamma, in_Z):
    func = defaultdict(int)
    for Z in [1, 6, 7, 8, 16]:
        with open("./data/Rayleigh/F-x-{}.dat".format(Z), "r") as f:
            data = f.readlines()[0:]
            x, y = zip(*[line.strip().split() for line in data])
            x = list(map(float, x))
            y = list(map(float, y))
            # 创建插值函数
            func[Z] = interpolate.interp1d(x, y)
    return Rayleigh_sigma_e(E_gamma, in_Z, func)


def draw_Rayleigh_curve():
    """
    无输入

    绘制液闪分子和各元素组分瑞利散射截面与入射光子能量的关系图

    返回(单位体积)(液闪分子康普顿效应截面, 液闪分子光电效应截面, Gamma能量取样点(begin, end, step))
    """

    for Z in [1, 6, 7, 8, 16]:
        with open("./data/Rayleigh/F-x-{}.dat".format(Z), "r") as f:
            data = f.readlines()[0:]
            x, y = zip(*[line.strip().split() for line in data])
            x = list(map(float, x))
            y = list(map(float, y))
            # 创建插值函数
            F[Z] = interpolate.interp1d(x, y)

    filename = "./data/LS.gdml"
    density_dict = extract_density(filename)
    Z_dict = Z_number(filename)

    # 给出一个gamma的能量取值列表
    E_gamma = np.linspace(0.00001, 10, 19999)
    plt.figure(figsize=(8, 10))
    plt.subplot(211)
    # 初始化液闪分子散射截面列表
    sigma_molecule = np.zeros_like(E_gamma)
    for element in Z_dict:
        sigma = []
        for e in E_gamma:
            Z = int(Z_dict[element])
            sigma.append(Rayleigh_sigma_e(e, Z, F))

        # 将当前图形添加到对应的子图中并计入液闪分子散射截面
        sigma_molecule += (
            np.array(sigma) * density_dict[element] / sum(density_dict.values())
        )
        plt.plot(E_gamma, sigma, label=f"{element}'s Rayleigh")
    Rayleigh_LSsigma = sigma_molecule
    plt.plot(E_gamma, sigma_molecule, label="LS's Rayleigh", linestyle="--", color="k")
    plt.xlabel("E_gamma")
    plt.ylabel("sigma")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.grid(True)
    plt.savefig("Rayleigh.pdf")


def integrate(func, x_start, x_end):
    """
    从 x_start 到 x_end, 积分f(x)^2
    """
    num_of_point = max(int((x_end - x_start) * 10), 200)
    if num_of_point == 0:
        num_of_point = 10
    x_val = np.linspace(x_start, x_end, num_of_point)
    f_x_val = func(x_val) ** 2
    return np.sum(f_x_val) * (x_end - x_start) / num_of_point


def calculate_b(i, kexi, PDF, X):
    return 1 - ((kexi[i + 1] - kexi[i]) / (X[i + 1] - X[i])) ** 2 / (
        PDF[i] * PDF[i + 1]
    )


def calculate_a(i, kexi, PDF, b, X):
    return ((kexi[i + 1] - kexi[i]) / (X[i + 1] - X[i])) / (PDF[i]) - b[i] - 1


def calculate_g(cos_theta):
    return (1 + cos_theta**2) / 2


def generate(X, A_coef, B_coef, kexi):
    """
    生成一个遵循X~kexi分布的随机数
    """
    interkexi = np.random.random()
    index = np.searchsorted(kexi, interkexi)
    if interkexi == kexi[index]:
        return X[index]

    kexi_pre, kexi_suf = kexi[index - 1], kexi[index]
    X_pre, X_suf = X[index - 1], X[index]
    A, B = A_coef[index - 1], B_coef[index - 1]

    nve = interkexi - kexi_pre
    delta = kexi_suf - kexi_pre

    ret = X_pre + (1 + A + B) * delta * nve / (
        delta**2 + A * delta * nve + B * (nve**2)
    ) * (X_suf - X_pre)

    # print(interkexi, X_pre, X_suf, ret, kexi_pre, kexi_suf)

    return ret


def sampleing_Rayleigh(E_gamma, Z):
    """
    传入光子能量(单位 MeV ), 相互作用质子的核电荷数

    返回光子散射角 theta 的 cos 形式, 也即 cos(theta)
    """
    with open("./data/Rayleigh/F-x-{}.dat".format(Z), "r") as f:
        data = f.readlines()[0:]
        x, y = zip(*[line.strip().split() for line in data])
        x = list(map(float, x))
        y = list(map(float, y))
        # 创建插值函数
        func = interpolate.interp1d(x, y)

    kappa = E_gamma / 0.511

    x_max = (20.6074 * 2 * kappa) ** 2
    x_min = 0

    x_test = np.linspace(x_min, x_max, 512)

    PDF_x_test = func(x_test) ** 2
    kexi_x_test = np.array([integrate(func, x_min, x_max_tmp) for x_max_tmp in x_test])
    total_integrate = integrate(func, x_min, x_max)

    kexi_norm = kexi_x_test / total_integrate  # 分布函数进行归一化处理
    PDF_x_norm = PDF_x_test / total_integrate  # 概率密度也归一化处理

    B_coef = np.array(
        [
            calculate_b(t, kexi_norm, PDF_x_norm, x_test)
            for t in range(x_test.shape[0] - 1)
        ]
    )  # 计算RITA算法中的B系数

    A_coef = np.array(
        [
            calculate_a(t, kexi_norm, PDF_x_norm, B_coef, x_test)
            for t in range(x_test.shape[0] - 1)
        ]
    )  # 计算RITA算法中的A系数

    while True:
        x_square = generate(x_test, A_coef, B_coef, kexi_norm)

        cos_theta = 1 - 1 / 2 * x_square / ((20.6074 * kappa) ** 2)

        kexi_test = np.random.random()

        if kexi_test <= calculate_g(cos_theta):
            return cos_theta


if __name__ == "__main__":
    draw_Rayleigh_curve()
    # print(sampleing(0.001, 1))
