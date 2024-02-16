.PHONY: all

# Gamma 光子能量, 单位 MeV
energy:=1
# 液闪物质中相互作用点的模拟次数
n_simulation:=2000 

all: pe_theta_$(energy) ce_photo_$(energy) cs_curve.pdf Rayleigh.pdf deposit.h5 deposit.html deposit_Bethe.h5 deposit_Bethe.html deposit_all_extend.h5 deposit_all_extend.html

# 光电效应电子出射角分布
pe_theta_$(energy):
	python3 photoelectric/photoelectric.py -e $(energy) -o $@


# 康普顿效应散射光子能量分布和散射角分布
ce_photo_$(energy):
	python3 compton/compton.py -e $(energy) -o $@


# 截面曲线
cs_curve.pdf:
	python3 ./interaction.py

# 瑞利散射截面曲线
Rayleigh.pdf:
	python3 ./Rayleigh/Rayleigh.py

# 模拟结果保存
deposit.h5:
	python3 ./main.py -e $(energy) -n_simulation $(n_simulation) 


# 考虑Bethe公式的模拟结果保存
deposit_Bethe.h5:
	python3 ./main_Bethe.py -e $(energy) -n_simulation $(n_simulation) 

# 考虑所有扩展信息的模拟结果保存
deposit_all_extend.h5:
	python3 ./main_all_extend.py -e $(energy) -n_simulation $(n_simulation)

# 沉积能量点分布
deposit.html deposit_dissymmetry.png: deposit.h5
	python3 ./analysis.py -i deposit.h5 

# 考虑Bethe公式的沉积能量点分布
deposit_Bethe.html deposit_Bethe_dissymmetry.png: deposit_Bethe.h5
	python3 ./analysis.py -i deposit_Bethe.h5

# 考虑所有扩展信息的沉积能量点分布
deposit_all_extend.html deposit_all_extend_dissymmetry.png: deposit_all_extend.h5
	python3 ./analysis.py -i deposit_all_extend.h5

clean:
	rm -f pe_theta_$(energy).pdf ce_photo_$(energy).pdf cs_curve.pdf deposit.html deposit_dissymmetry.png deposit.h5 deposit_Bethe.html deposit_Bethe_dissymmetry.png deposit_Bethe.h5 
	rm -f deposit_all_extend.h5 deposit_all_extend.html deposit_all_extend_dissymmetry.png Rayleigh.pdf

.DELETE_ON_ERROR:
.SECONDARY:
