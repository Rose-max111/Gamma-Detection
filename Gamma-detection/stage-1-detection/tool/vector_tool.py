import numpy as np


def rotate_vector(vector, axis, angle):
    """
    输入一个vector(3 dimensions numpy array), axis(3 dimensions numpy array), angle(float max = 360 degree)

    返回满足右手螺旋定则, 绕axis旋转angle的向量(3 dimensions numpy array)
    """
    # 角度转为弧度
    angle_rad = np.radians(angle)

    # 旋转轴的单位向量
    axis_unit = axis / np.linalg.norm(axis)

    # 旋转矩阵
    rotation_matrix = np.array(
        [
            [
                np.cos(angle_rad) + axis_unit[0] ** 2 * (1 - np.cos(angle_rad)),
                axis_unit[0] * axis_unit[1] * (1 - np.cos(angle_rad))
                - axis_unit[2] * np.sin(angle_rad),
                axis_unit[0] * axis_unit[2] * (1 - np.cos(angle_rad))
                + axis_unit[1] * np.sin(angle_rad),
            ],
            [
                axis_unit[0] * axis_unit[1] * (1 - np.cos(angle_rad))
                + axis_unit[2] * np.sin(angle_rad),
                np.cos(angle_rad) + axis_unit[1] ** 2 * (1 - np.cos(angle_rad)),
                axis_unit[1] * axis_unit[2] * (1 - np.cos(angle_rad))
                - axis_unit[0] * np.sin(angle_rad),
            ],
            [
                axis_unit[0] * axis_unit[2] * (1 - np.cos(angle_rad))
                - axis_unit[1] * np.sin(angle_rad),
                axis_unit[1] * axis_unit[2] * (1 - np.cos(angle_rad))
                + axis_unit[0] * np.sin(angle_rad),
                np.cos(angle_rad) + axis_unit[2] ** 2 * (1 - np.cos(angle_rad)),
            ],
        ]
    )

    # 计算旋转后的向量
    rotated_vector = np.dot(rotation_matrix, vector)
    return rotated_vector


def plot_test(vector, rotated_vector):
    """
    绘出vector与rotated_vector在三维坐标系中的示意图
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # 原始向量
    ax.quiver(
        0, 0, 0, vector[0], vector[1], vector[2], color="b", label="Original Vector"
    )

    # 旋转后的向量
    ax.quiver(
        0,
        0,
        0,
        rotated_vector[0],
        rotated_vector[1],
        rotated_vector[2],
        color="r",
        label="Rotated Vector",
    )

    # 设置坐标轴范围
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])

    # 设置坐标轴标签
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    # 添加图例
    ax.legend()

    # 显示图形
    plt.show()
