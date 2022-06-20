#anaconda py36excel

#xyzファイルのインポート
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import pathio

with open("C:/Users/kozo/Documents/PYTHON/Least_squares/jupyter_notebook/ati_1_1_out-Cloud.xyz", "r") as inf:
    iary = inf.readlines()

N = len(iary)
xyz = np.zeros((N, 3))

for i, ary in enumerate(iary):
    ary = ary.replace("\n", "")
    data = ary.split(" ")
    u, m, p = float(data[0]), float(data[1]), float(data[2])
    xyz[i] = [u, m, p]
    
#座標ごとの配列に転置
x, y, z = xyz.T

#座標ごとの最大値、最小値
ixyz = np.array([np.max(x), np.min(x), np.max(y), np.min(y), np.max(z), np.min(z)])
#座標ごとの長さ
lexyz = np.array([np.abs(ixyz[0] - ixyz[1]), np.abs(ixyz[2] - ixyz[3]), np.abs(ixyz[4] - ixyz[5])])

#yz平面データ
print("N: %d" %(N))
print("\npcd :")
print(xyz)

#グラフ描写
fig = plt.figure()
ax1 = Axes3D(fig)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('z')

ax1.scatter3D(x, y, z, label="Dataset", s=0.3)
plt.show()

#重心
G = np.zeros(3)
G = [np.average(x), np.average(y), np.average(z)]

print("\nG: (%f, %f, %f)" %(G[0], G[1], G[2]))

#各慣性モーメント
I = np.zeros(6) #ixx, iyy, izz, ixy, iyz, izx
for i in range(N):
    I[0] += (x[i] - G[0])**2
    I[1] += (y[i] - G[1])**2
    I[2] += (z[i] - G[2])**2
    I[3] += (x[i] - G[0])*(y[i] - G[1])
    I[4] += (y[i] - G[1])*(z[i] - G[2])
    I[5] += (z[i] - G[2])*(x[i] - G[0])

#慣性テンソル
A = np.zeros((3, 3))
A = [[I[0], I[3], I[5]],
     [I[3], I[1], I[4]],
     [I[5], I[4], I[2]]]

print("\n慣性テンソル:")
print(A[0])
print(A[1])
print(A[2])

#固有値Wと固有ベクトルV
W, V = np.zeros(3), np.zeros((3, 3))
W, V = np.linalg.eig(A)

print("\n固有値:")
print(W)
print("\n固有ベクトル:")
print(V)

#基準となるベクトル（固有値が最大となる固有ベクトル）
Widx = np.argmax(W)
V_std = np.array([V[Widx][0], V[Widx][1], V[Widx][2]])

#直線の方程式：x/a = y/b = z/c　を用いる
if np.argmax(lexyz) == 0:
    DN = lexyz[0]
    d = np.arange(int(ixyz[1]), int(ixyz[0]), 0.1)
    #print(d)
    st_pcd = np.zeros((len(d), 3))
    
    for i in range(len(d)):
        st_pcd[i] = [d[i], (d[i] - G[0])*V_std[1]/V_std[0] + G[1], (d[i] - G[0])*V_std[2]/V_std[0] + G[2]]
    
elif np.argmax(lexyz) == 1:
    print("同様")
elif np.argmax(lexyz) == 2:
    print("同様")

st_pcdt = st_pcd.T
#print(st_pcdt)

#グラフ描写
fig = plt.figure()
ax1 = Axes3D(fig)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('z')

ax1.scatter3D(x, y, z, label="Dataset", s=0.1)
ax1.scatter3D(st_pcdt[0], st_pcdt[1], st_pcdt[2], label="Dataset", s=5, color='red')
plt.legend()
plt.show()   


"""
ouf = pathio.outfle()
with open(ouf, "w") as f:
    np.savetxt(f, st_pcd)

"""
#輪切りにした時の各点群の重心Gnを求める

#分割数
n = 100
print("\n分割数: %d" %(n))

Gn = np.zeros((n, 3))

for i in range(n):
    xo, yo, zo =  st_pcd[int(len(d)*(i + 1)/n - len(d)/(2*n))]
    #ddは各平面までの最大距離の二乗
    dd = (st_pcd[int(len(d)*i/n)][0] - xo)**2 + (st_pcd[int(len(d)*i/n)][1] - yo)**2 + (st_pcd[int(len(d)*i/n)][2] - zo)**2
    nn = 0 #検出した点の数(初期化)
    
    #dd > (axo + byo + czo + d)^2 / a^2 + b^2 + c^2 を満たす点(xo, yo, zo)を検出する
    for j in range(N):
        ddi = (V_std[0]*x[j] + V_std[1]*y[j] + V_std[2]*z[j] - (V_std[0]*xo + V_std[1]*yo + V_std[2]*zo))**2 / (V_std[0]**2 + V_std[1]**2 + V_std[2]**2)
        
        if dd > ddi:
            Gn[i] += [x[j], y[j], z[j]]
            nn += 1
            
    Gn[i] /= nn

print("\n重心：")
print(Gn)

#グラフ表示
GnT = Gn.T
fig = plt.figure()
ax1 = Axes3D(fig)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('z')

ax1.scatter3D(x, y, z, label="Dataset", s=0.1)
ax1.scatter3D(st_pcdt[0], st_pcdt[1], st_pcdt[2], label="Dataset", s=5, color='red')
ax1.scatter3D(GnT[0], GnT[1], GnT[2], label="Dataset", s=50, color='black')
ax1.scatter3D(G[0], G[1], G[2], label="Dataset", s=50, color='yellow')
plt.legend()
plt.show()

ouf = pathio.outfle()
with open(ouf, "w") as f:
    np.savetxt(f, Gn)
