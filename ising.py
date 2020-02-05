import numpy as np
import matplotlib.pyplot as plt

def cal_delta_E(spin, i, j, Lattice, H, J):
    ## 计算翻转前后系统能量的变化
    top = i - 1
    bottom = (i + 1) % Lattice
    left = j - 1
    right = (j + 1) % Lattice
    E_loc_former = -1.0 * spin[i,j] * (J * (spin[top,j] + spin[bottom,j] + spin[i,left] + spin[i,right]) +  H )
    E_loc_next = spin[i,j] * (J * (spin[top,j] + spin[bottom,j] + spin[i,left] + spin[i,right]) +  H )
    delta_E = E_loc_next - E_loc_former
    return delta_E 
## 实验
sweeps = 1000 
relax = 50
K = 1
J = 1
H = 0 # 无外加磁场
results_dic = {} 
for Lattice in [5, 10, 20, 40, 80]: 
    # 初始化
    spin = np.random.random((Lattice,Lattice))
    spin[np.where(spin >= 0.5)] = 1
    spin[np.where(spin < 0.5)] = -1
    results = []
    for T in np.arange(0.1, 4.0, 0.05):
        mag = []
        # Monte Calro 每步翻转系统所有自旋
        for sweep in range(sweeps + relax):
            for i in range(Lattice):
                for j in range(Lattice):
                    delta_E = cal_delta_E(spin, i, j, Lattice, H, J)
                    if delta_E <= 0:
                        spin[i,j] = - spin[i,j]
                    # 几率翻转
                    elif np.exp((- delta_E)/(K * T)) > np.random.random():
                        spin[i,j] = - spin[i,j]
            if sweep % 100 == 0:
                print(sweep)
            M = abs(sum(sum(spin))) / (Lattice ** 2) # 平均磁矩正负等价，故取绝对值
            mag.append(M)
        result = sum(mag[relax:]) / sweeps # 统计平均值
        results.append(result)
    results_dic[Lattice] = results

## 作图
plt.figure(figsize=(12,9))
plt.subplot(231) 
plt.scatter(np.arange(0.1, 4.0, 0.05),results_dic[5],marker='o',c="",edgecolors='black',label='Lattice=5')
plt.grid(True)
plt.ylabel("M")
plt.legend()
plt.subplot(232) 
plt.scatter(np.arange(0.1, 4.0, 0.05),results_dic[10],marker='o',c="",edgecolors='blue',label='Lattice=10')
plt.grid(True)
plt.xlabel("T")
plt.legend()
plt.subplot(233) 
plt.scatter(np.arange(0.1, 4.0, 0.05),results_dic[20],marker='o',c="",edgecolors='red',label='Lattice=20')
plt.grid(True)
plt.legend()
plt.subplot(234) 
plt.scatter(np.arange(0.1, 4.0, 0.05),results_dic[40],marker='o',c="",edgecolors='green',label='Lattice=40')
plt.grid(True)
plt.legend()
plt.ylabel("M")
plt.subplot(235) 
plt.scatter(np.arange(0.1, 4.0, 0.05),results_dic[80],marker='o',c="",edgecolors='yellow',label='Lattice=80')
plt.grid(True)
plt.xlabel("T")
plt.legend()
plt.subplot(236) 
plt.scatter(np.arange(0.1, 4.0, 0.05),results_dic[5],marker='o',c="",edgecolors='black',label='Lattice=5')
plt.scatter(np.arange(0.1, 4.0, 0.05),results_dic[10],marker='o',c="",edgecolors='blue',label='Lattice=10')
plt.scatter(np.arange(0.1, 4.0, 0.05),results_dic[20],marker='o',c="",edgecolors='red',label='Lattice=20')
plt.scatter(np.arange(0.1, 4.0, 0.05),results_dic[40],marker='o',c="",edgecolors='green',label='Lattice=40')
plt.scatter(np.arange(0.1, 4.0, 0.05),results_dic[80],marker='o',c="",edgecolors='yellow',label='Lattice=80')
plt.grid(True)
plt.legend()
plt.savefig("T-Lattice.png")
plt.show()

##保存
f = open('results.txt','w')  
f.write(str(results_dic))  
f.close()
