#!/usr/bin/env python3
"""
    坐标系单位长度为mm
"""
import numpy as np
from tool.vector_tool import rotate_vector
from photoelectric.photoelectric import photoelectric
from compton.compton import compton_scattered_photon_sampling
from interaction import interacted
import h5py
import argparse
from scipy.constants import c, electron_mass, epsilon_0, elementary_charge
from interaction import extract_density, Z_number
import xml.etree.ElementTree as ET
import math
import random
from scipy.constants import c, elementary_charge

def Norm(Vec):
    total = np.sqrt(np.sum(Vec**2))
    return Vec / total

def A_number(filename):
    # 输入：文件名
    # 输出：一个字典。键值对为元素：摩尔质量A(g/mol)
    tree = ET.parse(filename)
    root = tree.getroot()
    string = extract_density(filename).keys()
    A_dict = {}

    for element in root.iter("element"):
        name = element.attrib.get("name")
        for prefix_length in range(1, len(name) + 1):
            prefix = name[:prefix_length]
            if prefix in string:
                break
        atom_value = float(element.find("atom").attrib.get("value"))
        A_dict[prefix] = atom_value

    return A_dict


def LS_density(filename):
    # 输入：文件名
    # 输出：一个浮点数（LS的密度g/cm3）
    tree = ET.parse(filename)
    root = tree.getroot()
    #定义XPath表达式以提取 <D> 元素中的 value 值
    d_xpath = ".//material[@name='LS']/D"

    # 提取 <D> 元素中的 value 值
    d_element = root.find(d_xpath)
    d_value = d_element.attrib["value"]

    # 输出 <D> 元素中的 value 值
    return(d_value)
def Bethe_x(E):
    """
    计算电子实际运动距离：

    输入: 电子能量(单个浮点数,单位:MeV)
    输出: 实际运动距离x(单个浮点数,单位:mm)
    """
    filename = "./data/LS.gdml"
    avogadro_number = 6.022e23
    LS_mass_density = float(LS_density(filename))

    density_dict = extract_density(filename)
    mol = sum(density_dict.values())
    for key in density_dict:
        density_dict[key] /= mol

    Z_dict = Z_number(filename)
    A_dict = A_number(filename)

    elec_energy = E * 1e6 * elementary_charge   # elec_energy单位：j
    step = 0.000001 # 单位：m
    e = elementary_charge

    # 使用Bethe公式计算电离损失
    def calculate_ion_loss(E, step, Z):

        # 计算洛伦兹因子γ
        gamma = E / (electron_mass * c**2) + 1

        # 通过洛伦兹因子γ计算速度v
        v = math.sqrt(1 - (1 / gamma**2)) * c
        beta = v / c
        loss_ion = (1 / (4 * np.pi * epsilon_0))**2 * ((2 * np.pi * e**4 * N* Z) / (electron_mass * v**2)) * (np.log((electron_mass * v**2 * E) / (2 * I**2 * (1 - beta**2))) - (np.log(2)*(2 * np.sqrt(1 - beta**2) - 1 + beta**2) + (1 - beta**2) + 1/8 * (1 - np.sqrt(1 - beta**2))**2))
        return loss_ion * step

    # 计算辐射损失
    def calculate_rad_loss(E, step, Z):
        loss_rad = ((1 / (4 * np.pi * epsilon_0))**2 * ((N * E * Z * (Z + 1) * e**4) / (137 * electron_mass * c**4)) *
                    (4 * np.log(2 * E / (electron_mass * c**2)) - 4 / 3))
        return loss_rad * step
    elements = list(density_dict.keys())
    probabilities = list(density_dict.values())
    Sum = 0
    N = 0
    for i in A_dict:
        Sum += A_dict[i] * density_dict[i]
    for i in density_dict:
        N += LS_mass_density * (density_dict[i] / Sum ) * avogadro_number * 1e6
    # 开始模拟
    x = 0
    while True:
        chosen_element = random.choices(elements, probabilities)[0]
        Z = Z_dict[chosen_element]
        I = Z * 10 * e
        elec_energy -= (calculate_rad_loss(elec_energy, step, Z) + calculate_ion_loss(elec_energy, step, Z))
        x += step
        if elec_energy <= 1e-30:
             return x * 1000


class main_atom:
    def __init__(self, Direction, Energy, Position):
        """
        输入 方向(3 Dimensions numpy array), 能量(float), 目前的位置(3 Dimensions numpy array)
        """
        self.dir = Norm(Direction)
        self.energy = Energy
        self.position = Position

    def next_atom(self, theta, phi, distance, new_energy):
        """
        输入 散射角theta(float), phi(float) —— 表示下一次移动相对当前Z轴的theta,phi角,
            移动长度distance(float),
            新粒子能量new_energy(float)

        返回一个main_atom类, 表示更新后的粒子各信息
        """
        y_axis_in_self = np.cross(np.array([0, 0, 1]), self.dir)

        vector_theta_fixed = rotate_vector(self.dir, y_axis_in_self, theta)
        vector_theta_phi_fixed = rotate_vector(vector_theta_fixed, self.dir, phi)

        vector_theta_phi_fixed = Norm(vector_theta_phi_fixed)

        return main_atom(
            vector_theta_phi_fixed,
            new_energy,
            self.position + vector_theta_phi_fixed * distance,
        )


if __name__ == "__main__":
    # 创建 ArgumentParser 对象
    parser = argparse.ArgumentParser()

    # 添加命令行参数
    parser.add_argument('-e', type=float)
    parser.add_argument('-n', type=int)
    parser.add_argument('-n_simulation', type=int)

    # 解析命令行参数
    args = parser.parse_args()

    # 获取参数值
    origin_E_gamma = args.e
    n_times = args.n_simulation# 随机次数
    origin_position = np.array([0, 0, 0])  # 初始反应位置
    origin_direction = np.array([1, 0, 0])  # 初始朝向
    atom_number = {1: 0, 6: 1, 7: 2, 8: 3, 16: 4}

    atom = [photoelectric(Z) for Z in [1, 6, 7, 8, 16]]

    output_data = []

    for T in range(n_times):
        act_photon = main_atom(origin_direction, origin_E_gamma, origin_position)

        '''
        测试：将光子初始方向随机 看看跟标准输出是否相似
        rand_phi = np.degrees(np.random.uniform(0, 2 * np.pi))
        rand_theta = np.degrees(np.random.uniform(0, np.pi))
        act_photon = act_photon.next_atom(rand_theta, rand_phi, 0, origin_E_gamma)
        '''


        Flag, Distance, act_atom_Z = interacted(act_photon.energy)
        output_data.append([])
        # output_data格式
        # 第一维度代表初始的第几个光子打入
        # 第二维度代表该光子衍生的第几次碰撞反应
        # 第三维度依次存储了 (初始是第几个光子(int)) (电子沉积坐标信息(3 dimensions)) (电子沉积能量(float))

        while True:
            # 把光子朝我们希望移动的地方移动
            act_photon = act_photon.next_atom(0, 0, Distance, act_photon.energy)

            if Flag == True:  # 进行光电效应
                cos_theta, elec_energy = atom[atom_number[act_atom_Z]].sampleing(
                    act_photon.energy
                )
                elec_theta = np.degrees(np.arccos(cos_theta))
                elec_phi = np.degrees(np.random.uniform(0, 2 * np.pi))  # 随机一个电子散射phi角

                if elec_energy <= 0:  # 如果模拟放出的电子能量 < 0 那么就是能量小于K层结合能，不放出电子
                    break

                end_elec = act_photon.next_atom(
                    elec_theta, elec_phi, Bethe_x(elec_energy), elec_energy
                )  # 计算出射电子的信息：在移动Bethe_x(elec_energy)后沉积、能量即为elec_energy
                output_data[T].append(
                    [
                        T,
                        end_elec.position[0],
                        end_elec.position[1],
                        end_elec.position[2],
                        end_elec.energy,
                    ]
                )
                break

            else:  # 进行康普顿效应
                epsilon, photon_theta, elec_theta = compton_scattered_photon_sampling(
                    act_photon.energy
                )
                photon_theta = np.degrees(photon_theta)  # 将photon返回的弧度值化为角度值
                elec_theta = np.degrees(elec_theta)  # 将electron返回的弧度值化为角度值

                elec_phi = np.degrees(np.random.uniform(0, 2 * np.pi))  # 随机一个电子散射phi角
                photon_phi = elec_phi + 180  # 光子在phi方向散射应与电子对称

                photon_energy = (
                    act_photon.energy * epsilon
                )  # 新光子能量等于入射光子能量 * (出射能量/入射能量 = epsilon)
                elec_energy = act_photon.energy - photon_energy  # 能量守恒
                distance = Bethe_x(elec_energy)
                end_elec = act_photon.next_atom(elec_theta, elec_phi, distance, elec_energy)
                output_data[T].append(
                    [
                        T,
                        end_elec.position[0],
                        end_elec.position[1],
                        end_elec.position[2],
                        end_elec.energy,
                    ]
                )  # 计算出射电子信息并计入输出

                next_photon = act_photon.next_atom(
                    photon_theta, photon_phi, 0, photon_energy
                )
                act_photon = next_photon  # 计算出射光子信息，调整出射光子的朝向(不移动)

            Flag, Distance, act_atom_Z = interacted(
                act_photon.energy
            )  # 采样下一次碰撞的类型(Flag)、下一次碰撞前移动的距离(Distance)、下一次碰撞的原子(act_atom_Z)
        print("fixed:", T)
    with h5py.File("deposit_Bethe.h5", "w") as opt:
        opt.create_group("/deposit_point")
        opt["deposit_point"].attrs["momentum"] = origin_E_gamma * elementary_charge * 1e6/ c
        opt["deposit_point"].attrs["direction"] = [1,0,0]
        for t in range(n_times):
            opt["deposit_point"][f"{t}"] = output_data[t]
