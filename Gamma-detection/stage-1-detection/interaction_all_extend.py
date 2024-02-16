from scipy.constants import c, electron_mass
from compton.compton import compton_sigma_e
from photoelectric.photoelectric import photoelectric
import math
import numpy as np
from matplotlib import pyplot as plt
import xml.etree.ElementTree as ET
import random
from Rayleigh.Rayleigh import get_Rayleigh_sigma

energy = electron_mass * c**2 / 1.6e-13  # 能量单位：MeV
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


def interacted(E_gamma):
    """
    输入光子能量E_gamma, 是否是第一次反应(默认为False)

    依次返回如下变量(共 3 个):

    Flag(Bool) 若Flag = False表示此次进行康普顿散射, 否则表示此次进行光电效应

    Distance(float)
    表示该光子到下一次反应位置处会运动的长度。若is_first = True则不需要考虑Distance造成的贝叶斯影响。单位为mm

    act_atom_Z(int) 表示该光子是与哪个原子发生的反应
    """
    filename = "./data/LS.gdml"
    density_dict = extract_density(filename)
    Z_dict = Z_number(filename)

    # compton
    compton_sigma_molecule = []
    # photoelectric
    photoelectric_sigma_molecule = []
    # Rayleigh
    Rayleigh_sigma_molecule = []

    for element in Z_dict:
        Z = int(Z_dict[element])

        photon = photoelectric(Z)

        compton_sigma = compton_sigma_e(E_gamma, Z)
        Rayleigh_sigma = get_Rayleigh_sigma(E_gamma, Z)
        photoelectric_sigma = photon.get_single_atom_photoelectric_cross_section(
            E_gamma
        )

        compton_sigma_molecule.append(compton_sigma * density_dict[element] / 1e24)
        photoelectric_sigma_molecule.append(
            photoelectric_sigma * density_dict[element] / 1e24
        )
        Rayleigh_sigma_molecule.append(Rayleigh_sigma * density_dict[element] / 1e24)

    sigma_molecule = (
        sum(compton_sigma_molecule)
        + sum(photoelectric_sigma_molecule)
        + sum(Rayleigh_sigma_molecule)
    )  # 单位体积内的总散射截面 (cm^-1)

    Lambda = 1 / sigma_molecule
    r = random.random()
    x = -Lambda * math.log(r) * 10  # 将抽样值转化为mm

    s = random.random()
    id = [6, 1, 8, 7, 16]  # 这里要根据输入文件的顺序 Z_dict中顺序是6 1 8 7 16

    prob_comp = sum(compton_sigma_molecule) / sigma_molecule  # 这是随机到康普顿散射的概率
    prob_photoelectric = (
        sum(photoelectric_sigma_molecule) / sigma_molecule
    )  # 这是随机到光电效应的概率
    prob_Rayleigh = sum(Rayleigh_sigma_molecule) / sigma_molecule  # 这是随机到瑞利散射的概率

    test_compton = compton_sigma_molecule / sum(
        compton_sigma_molecule
    )  # 此后为了抽样康普顿中哪种元素，需要对康普顿截面归一
    test_photoelectric = photoelectric_sigma_molecule / sum(
        photoelectric_sigma_molecule
    )  # 同理 也需要对光电截面归一
    test_Rayleigh = Rayleigh_sigma_molecule / sum(
        Rayleigh_sigma_molecule
    )  # 同理 也需要对瑞利截面归一

    if s <= prob_comp:  # 发生康普顿散射
        arr = np.random.multinomial(1, test_compton)
        for pos in range(5):
            if arr[pos] == 1:
                return 0, x, id[pos]
    elif s <= prob_comp + prob_photoelectric:  # 发生光电效应
        arr = np.random.multinomial(1, test_photoelectric)
        for pos in range(5):
            if arr[pos] == 1:
                return 1, x, id[pos]
    else:  # 发生瑞利散射
        arr = np.random.multinomial(1, test_Rayleigh)
        for pos in range(5):
            if arr[pos] == 1:
                return 2, x, id[pos]


def draw_curve():
    """
    无输入

    绘制液闪分子和各元素组分散射截面与入射光子能量的关系图

    返回(单位体积)(液闪分子康普顿效应截面, 液闪分子光电效应截面, Gamma能量取样点(begin, end, step))
    """
    filename = "./data/LS.gdml"
    density_dict = extract_density(filename)
    Z_dict = Z_number(filename)

    # 给出一个gamma的能量取值列表
    gamma_point = (0.001, 10, 9999)
    E_gamma = np.linspace(0.001, 10, 9999)
    plt.figure(figsize=(8, 10))
    plt.subplot(211)
    # 初始化液闪分子散射截面列表
    sigma_molecule = np.zeros_like(E_gamma)
    for element in Z_dict:
        sigma = []
        for e in E_gamma:
            Z = int(Z_dict[element])
            sigma.append(compton_sigma_e(e, Z))

        # 将当前图形添加到对应的子图中并计入液闪分子散射截面
        sigma_molecule += (
            np.array(sigma) * density_dict[element] / sum(density_dict.values())
        )
        plt.plot(E_gamma, sigma, label=f"{element}'s Compton")
    compton_LSsigma = sigma_molecule
    plt.plot(E_gamma, sigma_molecule, label="LS's Compton", linestyle="--", color="k")
    plt.xlabel("E_gamma")
    plt.ylabel("sigma")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.grid(True)
    plt.subplot(212)
    # 初始化液闪分子散射截面列表
    sigma_molecule = np.zeros_like(E_gamma)
    for element in Z_dict:
        sigma = []
        Z = int(Z_dict[element])
        photon = photoelectric(Z)
        for e in E_gamma:
            sigma.append(photon.get_single_atom_photoelectric_cross_section(e))

        # 将当前图形添加到对应的子图中并计入液闪分子散射截面
        sigma_molecule += (
            np.array(sigma) * density_dict[element] / sum(density_dict.values())
        )
        plt.plot(E_gamma, sigma, label=f"{element}'s Photoelectric")
    photoelectric_LSsigma = sigma_molecule
    plt.plot(
        E_gamma, sigma_molecule, label="LS's Photoelectric", linestyle="--", color="k"
    )
    plt.xlabel("E_gamma")
    plt.ylabel("sigma")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.grid(True)
    plt.subplots_adjust(top=0.9, bottom=0.1)
    plt.savefig("cs_curve.pdf", format="pdf")
    return (compton_LSsigma, photoelectric_LSsigma, gamma_point)


if __name__ == "__main__":
    draw_curve()
