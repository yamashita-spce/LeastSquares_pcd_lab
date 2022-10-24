# python=3.6.5 2022/08 yamashita

from distutils.sysconfig import get_config_var
import numpy as np
import os
import datetime
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import tkinter as tk
import tkinter.filedialog as fd
from scipy.optimize import minimize
from tqdm import tqdm

# グルーピングクラス
class PcdGroup:

    # gI ウェブ要素点群 gII 上フランジ要素点群 gIII 下フランジ要素点群
    gI = np.zeros((0,2))
    gII = np.zeros((0,2))
    gIII = np.zeros((0,2))

    # gp1: 上フランジ上面　gp2: 上フランジ下面　
    # gp3: ウェブ左面　gp4: ウェブ右面　
    # gp5: 下フランジ上面　gp6: 下フランジ下面
    gp1 = np.zeros((0,2))
    gp2 = np.zeros((0,2))
    gp3 = np.zeros((0,2))
    gp4 = np.zeros((0,2))
    gp5 = np.zeros((0,2))
    gp6 = np.zeros((0,2))

    # 各々の線分との距離
    gp1d = np.zeros(0)
    gp2d = np.zeros(0)
    gp3d = np.zeros(0)
    gp4d = np.zeros(0)
    gp5d = np.zeros(0)
    gp6d = np.zeros(0)

    def add_point(self, gp, p):

        return(np.append(gp, [p], axis=0))

    def gp_quart_range(self):

        self.gp1 = quart_range(self.gp1, self.gp1d)
        self.gp2 = quart_range(self.gp2, self.gp2d)
        self.gp3 = quart_range(self.gp3, self.gp3d)
        self.gp4 = quart_range(self.gp4, self.gp4d)
        self.gp5 = quart_range(self.gp5, self.gp5d)
        self.gp6 = quart_range(self.gp6, self.gp6d)
    
    def all_append(self):
        
        all_gp = np.zeros((0, 2))
        all_gp = np.append(all_gp, self.gp1, axis=0)
        all_gp = np.append(all_gp, self.gp2, axis=0)
        all_gp = np.append(all_gp, self.gp3, axis=0)
        all_gp = np.append(all_gp, self.gp4, axis=0)
        all_gp = np.append(all_gp, self.gp5, axis=0)
        all_gp = np.append(all_gp, self.gp6, axis=0)

        return all_gp
    
    def pcd_plot(self):

        print("\npcd plot...")
        plt.scatter(self.gp1.T[0], self.gp1.T[1], color="red")
        plt.scatter(self.gp2.T[0], self.gp2.T[1], color="blue")
        plt.scatter(self.gp3.T[0], self.gp3.T[1], color="green")
        plt.scatter(self.gp4.T[0], self.gp4.T[1], color="grey")
        plt.scatter(self.gp5.T[0], self.gp5.T[1], color="pink")
        plt.scatter(self.gp6.T[0], self.gp6.T[1], color="purple")
        plt.show()


# 主軸推定関数
def principal_axis_estimation(x, y, z):
    
    #重心
    G = np.array([np.average(x), np.average(y), np.average(z)])
    
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

    #固有値Wと固有ベクトルV
    W, V = np.zeros(3), np.zeros((3, 3))
    W, V = np.linalg.eig(A)
    V = V.T
    # V_inv = np.linalg.inv(V) #逆行列
    
    #基準となるベクトル（固有値が最大となる固有ベクトル）
    Widx = np.argmax(W)
    V_std = np.array([V[Widx][0], V[Widx][1], V[Widx][2]])
    # cp_V = V

    #固有値の小さい順のインデックス配列を取得
    w = np.argsort(W) 

    #法線ベクトルを直線に拡張（点群として） 直線の方程式：x/a = y/b = z/c　を用いる

    arasa = lexyz[np.argmax(lexyz)] / N #中心線を構成する点群粗さを定義
    
    if np.argmax(lexyz) == 0:
        d = np.arange(ixyz[1], ixyz[0], arasa) #主軸点群のデータ数
        st_pcd = np.zeros((len(d), 3)) # 主軸を構成する点群データ

        for i in range(len(d)):
            st_pcd[i] = [d[i], (d[i] - G[0])*V_std[1]/V_std[0] + G[1], (d[i] - G[0])*V_std[2]/V_std[0] + G[2]]

    elif np.argmax(lexyz) == 1:
        d = np.arange(ixyz[3], ixyz[2], arasa)
        st_pcd = np.zeros((len(d), 3))

        for i in range(len(d)):
            st_pcd[i] = [(d[i] - G[1])*V_std[0]/V_std[1] + G[0], d[i], (d[i] - G[1])*V_std[2]/V_std[1] + G[2]]

    elif np.argmax(lexyz) == 2:
        d = np.arange(ixyz[5], ixyz[4], arasa)
        st_pcd = np.zeros((len(d), 3))

        for i in range(len(d)):
            st_pcd[i] = [(d[i] - G[2])*V_std[0]/V_std[2] + G[0], (d[i] - G[2])*V_std[1]/V_std[2] + G[1], d[i]]

    if DESCRIPTION == True:
        #グラフ表示
        print("\npcd 3Dplot...")
        fig = plt.figure()
        ax1 = Axes3D(fig)
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_zlabel('z')

        ax1.scatter3D(x, y, z, label="Dataset", s=0.1)
        ax1.scatter3D(st_pcd.T[0], st_pcd.T[1], st_pcd.T[2], label="Dataset", s=5, color='red')
        plt.show()

    return V, V_std, w, st_pcd, d


# 断面ごとの二次元点群データを作成
def cross_section_pcd(V, V_std, x, y, z, st_pcd, dz):

    Sp = np.zeros((0, 2)) #二次元断面点群(座標変換後)

    #何分割目のものか -> dz
    xyzo = st_pcd[int(len(st_pcd)*dz/Nu)] #この点を通る断面を考える
    Z = np.dot(V, [[xyzo[0]],[xyzo[1]],[xyzo[2]]])[2][0] #系におけるZ座標

    #拾ってくる点の範囲を決定（デフォルトは推定300点拾うような範囲)
    Ds = (lexyz[np.argmax(lexyz)] * (D / 2)) / N
    
    #断面近郊の点を集める
    for i in range(N):
    # Nddi = np.sqrt((V_std[0]*(x[i] - xo) + V_std[1]*(y[i] - yo) + V_std[2]*(z[i] - zo))**2) / np.sqrt((V_std[0]**2 + V_std[1]**2 + V_std[2]**2)) #平面との距離
        Nddi = np.sqrt((V_std[0]*(x[i] - xyzo[0]) + V_std[1]*(y[i] - xyzo[1]) + V_std[2]*(z[i] - xyzo[2]))**2) #平面との距離
        if Nddi < Ds: 

            XYZ = np.dot(V, [[x[i]], [y[i]],[z[i]]])  #座標変換
            XY = np.array([XYZ[0][0], XYZ[1][0]])
            Sp = np.append(Sp, [XY], axis=0)

    return Sp, Z


#多変量最適化のための目的関数
def obj_func(p, pxy):
    #パラメータ: p = [Gx, Gy, h, Θ, α, β] pxy = [[x, y], bo]
    
    xyo = pxy[0]
    bo = pxy[1]
    SS = 0
    
    q, r = (bo/2)*np.cos(p[4]), (bo/2)*np.sin(p[4])
    s, t = (bo/2)*np.cos(p[5]), (bo/2)*np.sin(p[5])
      
    E = np.array([p[0] - (p[2]/2)*np.sin(p[3]), p[1] - (p[2]/2)*np.cos(p[3])])
    F = np.array([p[0] + (p[2]/2)*np.sin(p[3]), p[1] + (p[2]/2)*np.cos(p[3])])

    A = np.array([F[0] - q, F[1] + r])
    B = np.array([F[0] + q, F[1] - r])
    C = np.array([E[0] - s, E[1] + t])
    D = np.array([E[0] + s, E[1] - t])
    
    #ベクトル線分I, II, III の終始点行列を作成
    I = np.array([[E, F],
                  [A, B],
                  [C, D]])
    
    #各点においてベクトル線分の距離の二乗を計算し、その中で最小の値をSSに足す
    for j in range(len(xyo)):
        xy = xyo[j]
        ss = np.zeros(3)
        
        for i in range(3):
            
            #a, b は線分終始点配列
            a, b = I[i][0], I[i][1]
            d = b - a
            
            #線分のベクトル方程式:　a + td のtを求める
            t = np.dot((xy - a), d) / (np.linalg.norm(d)**2)

            #tに応じて、距離の二乗を求める
            if(t <= 0):
                ss[i] = np.linalg.norm(xy - a)**2
            elif(0 < t and t < 1):
                ss[i] = np.linalg.norm(xy - (a + d*t))**2
            elif(1 <= t):
                ss[i] = np.linalg.norm(xy - b)**2
            
        SS += np.min(ss)
       
    return SS

#結果と断面点群のプロット
def sectional_result_plot(Spp, P, bo):

    q, r = (bo/2)*np.cos(P.x[4]), (bo/2)*np.sin(P.x[4])
    s, t = (bo/2)*np.cos(P.x[5]), (bo/2)*np.sin(P.x[5])

    E = np.array([P.x[0] - (P.x[2]/2)*np.sin(P.x[3]), P.x[1] - (P.x[2]/2)*np.cos(P.x[3])])
    F = np.array([P.x[0] + (P.x[2]/2)*np.sin(P.x[3]), P.x[1] + (P.x[2]/2)*np.cos(P.x[3])])
    A = np.array([F[0] - q, F[1] + r])
    B = np.array([F[0] + q, F[1] - r])
    C = np.array([E[0] - s, E[1] + t])
    D = np.array([E[0] + s, E[1] - t])

    if DESCRIPTION == True:
        komakasa = 1000

        #ベクトル方程式 -->  [x, y] = a + t*(b - a)
        Iy = np.arange(E[1], F[1], (F[1] - E[1])/komakasa)
        IIx = np.arange(A[0], B[0],(B[0] - A[0])/komakasa)
        IIIx = np.arange(C[0], D[0], (D[0] - C[0])/komakasa)

        leI, leII, leIII = len(Iy), len(IIx), len(IIIx)
        I, II, III = np.zeros((leI, 2)), np.zeros((leII, 2)), np.zeros((leIII, 2))

        for i in range(leI):
            I[i][0] = E[0] + (i / leI)*(F[0] - E[0])
            I[i][1] = Iy[i]
        for i in range(leII):
            II[i][0] = IIx[i]
            II[i][1] = A[1] + (i / leII)*(B[1] - A[1])
        for i in range(leIII):
            III[i][0] = IIIx[i]
            III[i][1] = C[1] +(i / leIII)*(D[1] - C[1])    

        
        print("\npcd plot...")
        # plt.xlim(-1200, -1000) # x軸の表示範囲
        # plt.ylim(-780, -720)  # y軸の表示範囲
        plt.scatter(Spp.T[0], Spp.T[1], s=10)
        plt.scatter(I.T[0], I.T[1], s=1, color="red")
        plt.scatter(II.T[0], II.T[1], s=1, color="green")
        plt.scatter(III.T[0], III.T[1], s=1, color="orange")
        plt.scatter(P.x[0], P.x[1])
        plt.show()

    return(A, B, C, D, E, F)

# 推定重心から新たな座標系ベクトルの作成
def principal_axis_decision(G, Z, Q, w):

    L = np.zeros((Nu, 3))
    nV = np.zeros(((Nu, 3, 3))) #座標系ベクトル配列
    Ge = np.zeros((Nu, 3)) #重心配列
    Q_inv = np.linalg.inv(Q) #逆変換行列

    # 重心を元の座標系に戻す
    for i,g in enumerate(G):
        if np.all(g == 0):
            continue
        Ge[i] = np.dot(Q_inv, [g[0], g[1], Z[i]])

    # 主軸決定
    for i in range(Nu):

        if(i != 0 and i != (Nu - 1)):
            if(np.all(Ge[i] == 0)):
                continue
            elif(np.all(Ge[i - 1] == 0)):
                L[i] = (Ge[i + 1] - Ge[i]) / np.linalg.norm(Ge[i + 1] - Ge[i])
            elif(np.all(Ge[i + 1] == 0)):
                L[i] = (Ge[i] - Ge[i - 1]) / np.linalg.norm(Ge[i] - Ge[i - 1])
            else:
                L[i] = (Ge[i + 1] - Ge[i - 1]) / np.linalg.norm(Ge[i + 1] - Ge[i - 1])
        elif(i == 0):
            if(np.all(Ge[0]) == 0):
                continue
            else:
                L[0] = (Ge[1] - Ge[0]) / np.linalg.norm(Ge[1] - Ge[0])
        elif(i == (Nu - 1)):
            if(np.all(Ge[Nu -1]) == 0):
                continue
            else:
                L[Nu - 1] = (Ge[Nu - 1] - Ge[Nu - 2]) / np.linalg.norm(Ge[Nu - 1] - Ge[Nu - 2])

    # 座標系の決定
    for i in range(Nu):
        if(np.all(L[i]) == 0):
            continue
        else:
            # 右手座標系、かつ x: フランジ長さ方向, y: ウェブ長さ方向, z(奥行): 主軸方向になるように
            nV[i][0] = np.cross(L[i], Q[1])
            nV[i][1] = np.cross(nV[i][0], L[i])
            nV[i][2] = L[i]
    
    return nV

# グルーピング
def grouping(Sp, A):

    g = PcdGroup()

    #ベクトル線分I, II, III の終始点行列を作成
    I = np.array([[A[4], A[5]],
                [A[0], A[1]],
                [A[2], A[3]]])

    #ベクトル線分
    Ic = np.array([A[5] - A[4], A[1] - A[0], A[3] - A[2]])

    #各点においてベクトル線分の距離の二乗を計算
    for i in range(len(Sp)):
        xy = Sp[i]
        ss = np.zeros(3)
            
        for i in range(3):
                
            #a, b は線分終始点配列
            a, b = I[i][0], I[i][1]
            d = b - a
                
            #線分のベクトル方程式:　a + td のtを求める
            t = np.dot((xy - a), d) / (np.linalg.norm(d)**2)

            #tに応じて、距離の二乗を求める
            if(t <= 0):
                ss[i] = np.linalg.norm(xy - a)**2
            elif(0 < t and t < 1):
                ss[i] = np.linalg.norm(xy - (a + d*t))**2
            elif(1 <= t):
                ss[i] = np.linalg.norm(xy - b)**2
        
        #グルーピング
        if np.argmin(ss) == 0:
            # g.gI = np.append(g.gI, [xy], axis=0)
            g.gI = g.add_point(g.gI, xy)
        if np.argmin(ss) == 1:
            g.gII = np.append(g.gII, [xy], axis=0)
        if np.argmin(ss) == 2:
            g.gIII = np.append(g.gIII, [xy], axis=0)
    

    #さらに細かくグルーピング(gp)
    def distance_cal(a, b, xy):
        d = b - a
        # 線分のベクトル方程式:　a + td のtを求める
        t = np.dot((xy - a), d) / (np.linalg.norm(d)**2)
        
        return np.linalg.norm(xy - (a + d*t))**2
 
    for i in range(len(g.gI)):
        
        gd = distance_cal(I[0][0],  I[0][1], g.gI[i])
        
        if np.cross(Ic[0], g.gI[i] - A[4]) >= 0:
            g.gp3 = np.append(g.gp3, [g.gI[i]], axis=0)
            g.gp3d = np.append(g.gp3d, [gd])
        else:
            g.gp4 = np.append(g.gp4, [g.gI[i]], axis=0)
            g.gp4d = np.append(g.gp4d, [gd])
            
    for i in range(len(g.gII)):
        
        gd = distance_cal(I[1][0],  I[1][1], g.gII[i])
        
        if np.cross(Ic[1],g.gII[i] - A[0]) >= 0:
            g.gp1 = np.append(g.gp1, [g.gII[i]], axis=0)
            g.gp1d = np.append(g.gp1d, [gd])
        else:
            g.gp2 = np.append(g.gp2, [g.gII[i]], axis=0)  
            g.gp2d = np.append(g.gp2d, [gd])
            
    for i in range(len(g.gIII)):
        
        gd = distance_cal(I[2][0],  I[2][1], g.gIII[i])
        
        if np.cross(Ic[2], g.gIII[i] - A[2]) >= 0:
            g.gp5 = np.append(g.gp5, [g.gIII[i]], axis=0)
            g.gp5d = np.append(g.gp5d, [gd])
        else:
            g.gp6 = np.append(g.gp6, [g.gIII[i]], axis=0) 
            g.gp6d = np.append(g.gp6d, [gd])
    
    if DESCRIPTION == True: g.pcd_plot()

    return g

# 四分位範囲の点群を抜き出す関数
def quart_range(p, d):
    sp =  np.zeros((0, 2))
    sidx = np.argsort(d)

    if len(p) >= 4:
        for j in sidx[int(len(sidx)*1/4):int(len(sidx)*3/4)]:
            sp = np.append(sp, [p[j]], axis=0)
        
    return sp

# フランジ線分を決定するための目的関数
def obj_func_flange(ps, xyxy):

    # y = ax + b, y = ax + c　において
    #パラメータps = [a, b, c], xyxy = [[xy],[xy]]の三次元配列
    xy1 = xyxy[0]
    xy2 = xyxy[1]
    Ss = 0

    for i in range(len(xy1)):
        Ss += (ps[0]*xy1[i][0] + ps[1] - xy1[i][1])**2
    
    for i in range(len(xy2)):
        Ss += (ps[0]*xy2[i][0] + ps[2] - xy2[i][1])**2
        
    return Ss

# ウェブ線分を決定するための目的関数
def obj_func_web(ps, xyxy):

    # x = ay + b, x = ay + c　において
    #パラメータps = [a, b, c], xyxy = [[xy],[xy]]の三次元配列
    xy1 = xyxy[0]
    xy2 = xyxy[1]
    Ss = 0

    for i in range(len(xy1)):
        Ss += (ps[0]*xy1[i][1] + ps[1] - xy1[i][0])**2
    
    for i in range(len(xy2)):
        Ss += (ps[0]*xy2[i][1] + ps[2] - xy2[i][0])**2
        
    return Ss

def uniform_flange_thickness(pII, pIII):
    
    #フランジの厚さの平均
    L = (np.abs(pII[1] - pII[2])*np.cos(np.arctan(pII[0])) + np.abs(pIII[1] - pIII[2])*np.cos(np.arctan(pIII[0])))/2   
    
    #フランジの厚さが平均になるように線分を平行移動
    dyII = (L / np.cos(np.arctan(pII[0]))) - np.abs(pII[1] - pII[2]) 
    dyIII = (L / np.cos(np.arctan(pIII[0]))) - np.abs(pIII[1] - pIII[2]) 
    pII[1] += dyII/2 
    pII[2] -= dyII/2
    pIII[1] += dyIII/2
    pIII[2] -= dyIII/2

    return pII, pIII, L

#直線からの距離が(2/5)L の範囲にあり、Iの中心線における点との最大距離×５/4以上の距離である点を抜き取りと右側の点群に分ける
def side_flange(P, Iy, Ic, a, b, L): 

    #y = ax + b, yはIのyの係数、IcはＩの平均ｘ切片
    Lp = np.zeros((0,2))
    Rp = np.zeros((0,2))
    Id = np.zeros(len(P))
    Idmax = 0
    
    for i,p in enumerate(P):      
        Id[i] = np.abs(p[0] - Iy*p[1] - Ic) / np.sqrt(Iy**2 + 1)
        
        if Id[i] > Idmax:
            Idmax = Id[i]
    
    for i,p in enumerate(P):
        IId = np.abs(-a*p[0] + p[1] - b) / np.sqrt(a**2 + 1)
        
        if IId <= L*(2 / 5) and Id[i] > Idmax*4/5:
            if p[1]*Iy+Ic-p[0] > 0:
                Rp = np.append(Rp, [p], axis=0)
            else:
                Lp = np.append(Lp, [p], axis=0)
    
    return Rp, Lp

# フランジ横線分を決定するための目的関数
def obj_func_side_flange(c, pc):

    #x = -IId*y + c[0], x = -IId*y + c[1]　（IIdはIIの傾き,もしくはIIIの傾き)
    #pc = [[Lp], [Rp], IId]を格納する配列 (Lp, Rpは二次元点群データ)
    Lp = pc[0]
    Rp = pc[1]
    Ss = 0
    
    for lp in Lp:
        Ss += (-pc[2]*lp[1] + c[0] - lp[0])**2
    for rp in Rp:
        Ss += (-pc[2]*rp[1] + c[1] - rp[0])**2
    
    return Ss

def uniform_side_flange(Pt, Pb, PII, PIII):

    #フランジ長さの平均
    U = np.abs(((Pt[1] - Pt[0])*np.cos(np.arctan(PII[0])) + (Pb[1] - Pb[0])*np.cos(np.arctan(PIII[0]))) / 2)
    
    #フランジの長さが平均になるように線分を平行移動
    dxII = (U / np.cos(np.arctan(PII[0]))) - np.abs(Pt[1] - Pt[0]) 
    dxIII = (U / np.cos(np.arctan(PIII[0]))) - np.abs(Pb[1] - Pb[0]) 
    Pt[0] -= dxII/2 
    Pt[1] += dxII/2
    Pb[0] -= dxIII/2
    Pb[1] += dxIII/2

    return Pt, Pb, U

# webを構成する頂点から重心点を推定
def cogravity_est(PI, PII, PIII):

    vtx_web = np.zeros((4, 2))
    vtx_web[0] = slv_equt(PIII[0], PIII[1], PI[0], PI[1]) #E_left
    vtx_web[1] = slv_equt(PIII[0], PIII[1], PI[0], PI[2]) #E_right
    vtx_web[2] = slv_equt(PII[0], PII[2], PI[0], PI[1]) #F_left
    vtx_web[3] = slv_equt(PII[0], PII[2], PI[0], PI[2]) #F_right

    mean = np.array([np.mean(vtx_web.T[0]), np.mean(vtx_web.T[1])]) #web座標平均

    return mean


# 頂点を決定する
def vertex_coordinate_detection(PI, PII, PIII, Pts, Pbs):
    
    vtx_p = np.zeros((12, 2))
    vtx_p[0] = slv_equt(PII[0], PII[1], -PII[0], Pts[0]) #A_up
    vtx_p[1] = slv_equt(PII[0], PII[2], -PII[0], Pts[0]) #B_up
    vtx_p[2] = slv_equt(PII[0], PII[1], -PII[0], Pts[1]) #A_down
    vtx_p[3] = slv_equt(PII[0], PII[2], -PII[0], Pts[1]) #B_down
    vtx_p[4] = slv_equt(PIII[0], PIII[1], -PIII[0], Pbs[0]) #C_up
    vtx_p[5] = slv_equt(PIII[0], PIII[2], -PIII[0], Pbs[0]) #D_up
    vtx_p[6] = slv_equt(PIII[0], PIII[1], -PIII[0], Pbs[1]) #C_down
    vtx_p[7] = slv_equt(PIII[0], PIII[2], -PIII[0], Pbs[1]) #D_down
    vtx_p[8] = slv_equt(PIII[0], PIII[1], PI[0], PI[1]) #E_left
    vtx_p[9] = slv_equt(PIII[0], PIII[1], PI[0], PI[2]) #E_right
    vtx_p[10] = slv_equt(PII[0], PII[2], PI[0], PI[1]) #F_left
    vtx_p[11] = slv_equt(PII[0], PII[2], PI[0], PI[2]) #F_right

    # 重心座標の取得
    G = np.array([np.mean(vtx_p.T[0][8:12]), np.mean(vtx_p.T[1][8:12])])

    return vtx_p, G

#y = [yd]*x + [yc], x = [xd]*y + [xc]の解を求める関数
def slv_equt(yd, yc, xd, xc): 
    y = (yc + yd*xc) / (1 - yd*xd)
    x = xd*y + xc  
    return x, y

# 輪郭を構成する点を作成
def create_contour_points(vtx_p):

    def line_segment(p1, p2): #p1, p2は始点と終点
        P = np.zeros((pN, p1.shape[0]))
        
        for i in range(len(P)):
            P[i] = p1 + ((i + 1)/pN)*(p2 - p1)
        
        return P

    seg_p = np.zeros((0, 2))
    seg_p = np.append(seg_p, line_segment(vtx_p[0], vtx_p[1]), axis=0)
    seg_p = np.append(seg_p, line_segment(vtx_p[2], vtx_p[3]), axis=0)
    seg_p = np.append(seg_p, line_segment(vtx_p[4], vtx_p[5]), axis=0)
    seg_p = np.append(seg_p, line_segment(vtx_p[6], vtx_p[7]), axis=0)
    seg_p = np.append(seg_p, line_segment(vtx_p[8], vtx_p[10]), axis=0)
    seg_p = np.append(seg_p, line_segment(vtx_p[9], vtx_p[11]), axis=0)
    seg_p = np.append(seg_p, line_segment(vtx_p[0], vtx_p[2]), axis=0)
    seg_p = np.append(seg_p, line_segment(vtx_p[1], vtx_p[10]), axis=0)
    seg_p = np.append(seg_p, line_segment(vtx_p[11], vtx_p[3]), axis=0)
    seg_p = np.append(seg_p, line_segment(vtx_p[4], vtx_p[8]), axis=0)
    seg_p = np.append(seg_p, line_segment(vtx_p[9], vtx_p[6]), axis=0)
    seg_p = np.append(seg_p, line_segment(vtx_p[5], vtx_p[7]), axis=0)

    return seg_p

# 座標系を元に戻す
def retransform_coordinate_system(p, V, Z):

    vp = np.zeros((len(p), 3))
    inv_V = np.linalg.inv(V)

    for i in range(len(p)):
        vp[i] = np.dot(inv_V, [[p[i][0]],[p[i][1]],[Z]]).T

    return vp

# 重心座標を元の座標系に戻す
def rc_gv(G, V, Z):
    inv_V = np.linalg.inv(V)

    return np.dot(inv_V, [[G[0]],[G[1]],[Z]]).T

# 分割数の指定
def division_number_input():

    global Nu
    while True:
        print("分割数の入力 > ", end="")
        try:
            Nu = int(input())
            break
        except:
            print("[ERROR] 整数を入力してください\n")
            continue
    
    if Nu >= N: 
        print("[EEROR] 分割数が多すぎます。'%d'以下の分割数を指定してください" %(N))
        exit()

# ログ出力
def output_log(txt, out):

    if out == True:
        print(txt)
    with open(epath, "a") as ef:
        ef.write(txt + "\n")

#インポートファイルのpathの取得
def infle():

    print("Choose input file ...")
    root = tk.Tk()
    root.withdraw()
    inf = fd.askopenfilename(
        filetypes=[("xyz", ".xyz"), ("txt", ".txt")],
    )
    print("Input file : %s\n" %(inf))

    return(inf)

#アウトプットファイルを選択・作成
def outfle():
    global ouf_vtx, ouf_cont, ouf_td

    print("Create output file ...")
    root = tk.Tk()
    root.withdraw()
    ouf = fd.asksaveasfilename(
        title = "名前を付けて保存",
        filetypes=[("xyz", ".xyz")],
        initialdir=dpath,
        defaultextension = "xyz"
    )

    ouf_vtx = os.path.splitext(ouf)[0] + "_vtx.xyz" #頂点座標を出力するファイル
    ouf_cont = os.path.splitext(ouf)[0] + "_cont.xyz" #輪郭線分を点群で出力するファイル
    ouf_td = os.path.splitext(ouf)[0] + "_twodimension.tmp" #二次元系の頂点座標を出力するファイル

    vf = open(ouf_vtx, "w")
    cf = open(ouf_cont, "w")
    tf = open(ouf_td, "w")
    vf.close()
    cf.close()
    tf.close()

    output_log("\nOutput file : %s" %(ouf_vtx), True)
    output_log("Output file : %s" %(ouf_cont), True)
    output_log("Output file : %s\n" %(ouf_td), True)

#重心点plot
def grav_plot(p, g):

    print("grav plot...")
    plt.scatter(p.T[0], p.T[1])
    plt.scatter(g[0], g[1], s=40, color="red", marker="*")
    plt.show()

# 二次元plot
def sur_plot(p, sp):
    
    print("plot...")
    plt.scatter(sp.T[0], sp.T[1])
    plt.scatter(p.T[0], p.T[1], s=0.5, color="red")
    plt.show()

# 3次元plot
def V_plot(p):
    #グラフ表示
    print("\npcd 3Dplot...")
    fig = plt.figure()
    ax1 = Axes3D(fig)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')

    ax1.scatter3D(p.T[0], p.T[1], p.T[2], label="Dataset")
    plt.show()

# パラメータのアウトプット、ログファイルへの書き込み
def parameter_output(nn, gx, gy, h, d, du, dd, bo, string, number):
    
    output_log("\n%s step = %d" %(string, number), DETAIL)
    output_log("==========================================", DETAIL)
    output_log("点群総数：%d" %(nn), DETAIL)
    output_log("重心点 : ( %.4f, %.4f )" %(gx, gy), DETAIL)
    output_log("高さ : %.4f" %(h), DETAIL)
    output_log("座標軸との傾き : %.4f[rad]" %(d), DETAIL)
    output_log("ウェブと上フランジの傾き : %.4fπ[rad]" %(0.5 + du/np.pi), DETAIL)
    output_log("ウェブと下フランジの傾き : %.4fπ[rad]" %(0.5 + du/np.pi), DETAIL)
    output_log("幅（固定値）: %.4f" %(bo), DETAIL)
    output_log("==========================================", DETAIL)

def line_segment_output(pI, pII, pIII, pts, pbs, number):

    output_log("\n線分推定結果 step = %d" %(number), DETAIL)
    output_log("\n============================", DETAIL)
    output_log("I': x = %.3fy + %.3f" %(pI[0], pI[1]), DETAIL)
    output_log("I'': x = %.3fy + %.3f" %(pI[0], pI[2]), DETAIL)
    output_log("II': y = %.3fx + %.3f" %(pII[0], pII[1]), DETAIL)
    output_log("II'': y = %.3fx + %.3f" %(pII[0], pII[2]), DETAIL)
    output_log("III': y = %.3fx + %.3f" %(pIII[0], pIII[1]), DETAIL)
    output_log("III'': y = %.3fx + %.3f" %(pIII[0], pIII[2]), DETAIL)
    output_log("IIL : x = %.3fy + %.3f" %(-pII[0], pts[0]), DETAIL)
    output_log("IIR : x = %.3fy + %.3f" %(-pII[0], pts[1]), DETAIL)
    output_log("IIIL : x = %.3fy + %.3f" %(-pIII[0], pbs[0]), DETAIL)
    output_log("IIIR : x = %.3fy + %.3f" %(-pIII[0], pbs[1]), DETAIL)
    output_log("==============================", DETAIL)


def main(): 

    global  DETAIL, DESCRIPTION, N, Nu, pN, D, input_file, ouf_vtx, ouf_cont, ouf_td, epath, dpath, ixyz, lexyz

    DETAIL = False #詳細を出力する
    DESCRIPTION = False #点群を描写する
    D = 300 #分割したときに拾う点群の数の指定（おおよそ）
    pN = 1000 #輪郭を構成する点の各線分における点群の数（輪郭点不要の場合は少なくしたほうが良い）

    try:
        # 点群ファイルの読込
        input_file = infle()
        with open(input_file, "r") as inf:
            iary = inf.readlines()
        N = len(iary)
        xyz = np.zeros((N, 3))

        for i, ary in enumerate(iary):
            ary = ary.replace("\n", "")
            data = ary.split(" ")
            u, m, p = float(data[0]), float(data[1]), float(data[2])
            xyz[i] = [u, m, p]

        # ログファイルの作成
        dpath, adder = os.path.split(input_file)
        dt_now = datetime.datetime.now()
        epath = dpath + "/" + str(dt_now.strftime('%Y%m%d-%H%M%S')) + "_log.txt"
        etxt = open(epath, mode='w')
        etxt.close()

        output_log("Input file : %s" %(input_file), False)
        output_log("Total number of pcd : %d\n" %(len(xyz)), True)
        output_log("log file path : %s\n" %(epath), True)

        # 書き込みファイルの作成、指定
        outfle()
        # 分割数の指定
        division_number_input()
        #座標ごとの配列に転置
        x, y, z = xyz.T
        #座標ごとの最大値、最小値
        ixyz = np.array([np.max(x), np.min(x), np.max(y), np.min(y), np.max(z), np.min(z)])
        #座標ごとの点群の最大距離
        lexyz = np.array([np.abs(ixyz[0] - ixyz[1]), np.abs(ixyz[2] - ixyz[3]), np.abs(ixyz[4] - ixyz[5])])
       
    except Exception as e:
        print(e)
        exit()


    # =========================================================主軸推定=========================================================
    
    V, V_std, w, st_pcd, d = principal_axis_estimation(x, y, z)
    Q = np.array([V[w[0]], V[w[1]], V[w[2]]]) #変換行列（x: フランジ長さ方向, y: ウェブ長さ方向, z: 主軸方向になるような変換行列にする)

    # 重心推定
    output_log("< center of gravity estimation >", False)
    ZZ, Gn = np.zeros(Nu), np.zeros((Nu, 2))
    pber1 = tqdm(total = Nu)

    for i in range(Nu):
        pber1.set_description("center of gravity estimation")

        # 断面点群の取得（二次元系）
        Sp, ZZ[i] = cross_section_pcd(Q, V_std, x, y, z, st_pcd, i)

        #断面点群が存在しない場合
        if len(Sp) == 0: 
            output_log("\n[WARNING] '%d' Point cloud does not exist near reference cross section " %(i), True)
            Gn[i][0], Gn[i][1] = 0, 0
            pber1.update(1)
            continue
        
        if DESCRIPTION == True: plt.scatter(Sp.T[0], Sp.T[1]); plt.show()

        # パラメータ初期値の設定　p = [Gx, Gy, h, rad] pxy = [x, y]
        p = np.array([np.average(Sp.T[0]), np.average(Sp.T[1]), np.abs(np.max(Sp.T[1]) - np.min(Sp.T[1])), 0, 0, 0])
        bo = np.abs(np.max(Sp.T[0]) - np.min(Sp.T[0]))
        parameter_output(len(Sp), p[0], p[1], p[2], 0, 0, 0, bo, "初期値", i)

        # 最適化
        P = minimize(obj_func, p, args=[Sp, bo], method="powell")
        parameter_output(len(Sp), P.x[0], P.x[1], P.x[2], P.x[3], P.x[4], P.x[5], bo, "最適化", i)

        A = sectional_result_plot(Sp, P, bo) #頂点の取得
        g = grouping(Sp, A) #グルーピング

        #　実験用
        # with open("data2_group_point_exam%d.txt" %(i), "w") as f:
        #     f.write("%d %d %d %d %d %d\n" %(len(g.gp1), len(g.gp2), len(g.gp3), len(g.gp4), len(g.gp5), len(g.gp6)))
        #     np.savetxt(f ,A)
        #     np.savetxt(f, g.gp1)
        #     np.savetxt(f, g.gp2)
        #     np.savetxt(f, g.gp3)
        #     np.savetxt(f, g.gp4)
        #     np.savetxt(f, g.gp5)
        #     np.savetxt(f, g.gp6)
        #     np.savetxt(f, g.gp1d)
        #     np.savetxt(f, g.gp2d)
        #     np.savetxt(f, g.gp3d)
        #     np.savetxt(f, g.gp4d)
        #     np.savetxt(f, g.gp5d)
        #     np.savetxt(f, g.gp6d)     

        g.gp_quart_range() #四分位範囲を抜き出す
        
        #上フランジパラメータ初期値
        ds = (A[1][1] - A[0][1])/(A[1][0] - A[0][0])
        ps = np.array([ds, A[0][1] - ds*A[0][0], A[0][1] - ds*A[0][0]])
        PII = minimize(obj_func_flange, ps, args=[g.gp1, g.gp2], method="powell") #最適化

        #ウェブパラメータ初期値
        ds = (A[5][0] - A[4][0])/(A[5][1] - A[4][1])
        ps = np.array([ds, A[4][0] - ds*A[4][1], A[4][0] - ds*A[4][1]])
        PI = minimize(obj_func_web, ps, args=[g.gp3, g.gp4], method="powell") #最適化

        #下フランジパラメータ初期値
        ds = (A[4][1] - A[3][1])/(A[4][0] - A[3][0])
        ps = np.array([ds, A[3][1] - ds*A[3][0], A[3][1] - ds*A[3][0]])
        PIII = minimize(obj_func_flange, ps, args=[g.gp5, g.gp6], method="powell") #最適化
        
        PII.x, PIII.x, L = uniform_flange_thickness(PII.x, PIII.x) #フランジ厚さを調整

        Gn[i][0], Gn[i][1] = cogravity_est(PI.x, PII.x, PIII.x) #重心点の推定  

        if DESCRIPTION == True: grav_plot(Sp, Gn[i]) #重心点の確認

        # Gn[i][0], Gn[i][1] = P.x[0], P.x[1]
        pber1.update(1)

    pber1.close()

    rV = principal_axis_decision(Gn, ZZ, Q, w) #各断面の主軸に対する二次元変換行列を作成
    

    # =========================================================フィッティング・特徴点取得=========================================================

    output_log("< Fitting / Feature point acquisition >", False)
    vtx, cont = np.zeros((0, 3)), np.zeros((0, 3))

    rZZ = np.zeros(Nu) #系におけるz座標
    pber2 = tqdm(total = Nu)

    for i in range(Nu):
        pber2.set_description("Fitting / Feature point acquisition")
        if(np.all(rV[i]) == 0):
            pber2.update(1)
            continue
        
        # 断面点群の取得（二次元系）
        rSp, rZZ[i] = cross_section_pcd(rV[i], rV[i][2], x, y, z, st_pcd, i)
        if DESCRIPTION == True: plt.scatter(rSp.T[0], rSp.T[1]); plt.show()

        # パラメータ初期値の設定　p = [Gx, Gy, h, rad] pxy = [x, y]
        p = np.array([np.average(rSp.T[0]), np.average(rSp.T[1]), np.abs(np.max(rSp.T[1]) - np.min(rSp.T[1])), 0, 0, 0])
        bo = np.abs(np.max(rSp.T[0]) - np.min(rSp.T[0]))
        parameter_output(len(rSp), p[0], p[1], p[2], 0, 0, 0, bo, "初期値", i)

        # 最適化
        P = minimize(obj_func, p, args=[rSp, bo], method="powell")
        parameter_output(len(rSp), P.x[0], P.x[1], P.x[2], P.x[3], P.x[4], P.x[5], bo, "フィッティング後パラメータ", i)
        
        A = sectional_result_plot(rSp, P, bo) #頂点の取得
        rg = grouping(rSp, A) #グルーピング
        rg.gp_quart_range() #四分位範囲を抜き出す
        if DESCRIPTION == True: rg.pcd_plot(); print(i)
        
        #上フランジパラメータ初期値
        ds = (A[1][1] - A[0][1])/(A[1][0] - A[0][0])
        ps = np.array([ds, A[0][1] - ds*A[0][0], A[0][1] - ds*A[0][0]])
        PsII = minimize(obj_func_flange, ps, args=[rg.gp1, rg.gp2], method="powell") #最適化
        
        #ウェブパラメータ初期値
        ds = (A[5][0] - A[4][0])/(A[5][1] - A[4][1])
        ps = np.array([ds, A[4][0] - ds*A[4][1], A[4][0] - ds*A[4][1]])
        PsI = minimize(obj_func_web, ps, args=[rg.gp3, rg.gp4], method="powell") #最適化
        
        #下フランジパラメータ初期値
        ds = (A[4][1] - A[3][1])/(A[4][0] - A[3][0])
        ps = np.array([ds, A[3][1] - ds*A[3][0], A[3][1] - ds*A[3][0]])
        PsIII = minimize(obj_func_flange, ps, args=[rg.gp5, rg.gp6], method="powell") #最適化
        
        PsII.x, PsIII.x, L = uniform_flange_thickness(PsII.x, PsIII.x) #フランジ厚さを調整
        
        #上フランジの横面を構成する点
        tpLp, tpRp = side_flange(rg.gII, PsI.x[0], (PsI.x[1] + PsI.x[2])/2, PsII.x[0], (PsII.x[1]+PsII.x[2])/2, L)
        #下フランジの横面を構成する点
        bmLp, bmRp = side_flange(rg.gIII, PsI.x[0], (PsI.x[1] + PsI.x[2])/2, PsIII.x[0], (PsIII.x[1]+PsIII.x[2])/2, L)
        
        if len(tpLp) == 0 or len(tpRp) == 0 or len(bmLp) == 0 or len(bmRp) == 0:
            output_log("[WARNING] '%d' Point cloud for the flange side does not exist. Resolved by increasing the number of point clouds" %(i), True) 
            pber2.update(1)
            continue

        #上フランジ横面パラメータ初期値
        tc = np.array([tpLp[0][0], tpRp[0][0]])
        Pts = minimize(obj_func_side_flange, tc, args=[tpLp, tpRp, PsII.x[0]], method="powell") #最適化   

        #下フランジ横面パラメータ初期値
        tc = np.array([bmLp[0][0], bmRp[0][0]])
        Pbs = minimize(obj_func_side_flange, tc, args=[bmLp, bmRp, PsIII.x[0]], method="powell") #最適化   

        Pts.x, Pbs.x, U = uniform_side_flange(Pts.x, Pbs.x, PsII.x, PsIII.x) #フランジ幅の調整
        line_segment_output(PsI.x, PsII.x, PsIII.x, Pts.x, Pbs.x, i)

        s_vtx, Gt = vertex_coordinate_detection(PsI.x, PsII.x, PsIII.x, Pts.x, Pbs.x) #頂点を決定  
        s_seg = create_contour_points(s_vtx) #輪郭点を作成
        if DESCRIPTION == True: sur_plot(s_seg, rSp)

        vtx = np.append(vtx, retransform_coordinate_system(s_vtx, rV[i], rZZ[i]), axis=0)
        cont = np.append(cont, retransform_coordinate_system(s_seg, rV[i], rZZ[i]), axis=0)
        if DESCRIPTION == True: V_plot(cont)

        # dat_file生成用tmpファイルに書き込み
        with open(ouf_td, 'a') as f_handle:
            np.savetxt(f_handle, rc_gv(Gt, rV[i], rZZ[i])) #三次元重心座標の保存
            np.savetxt(f_handle, s_vtx) #データの保存
            

        pber2.update(1)

    pber2.close()
    
    # データの保存
    np.savetxt(ouf_vtx, vtx)
    np.savetxt(ouf_cont, cont)
    


if __name__ == "__main__":
    print("===============================================")
    print("noIncompletePcd.py   ( version: 1.1.1 ) ")
    print("python 3.6.5, anaconda 4.10.3 ")
    print("Lib: numpy:1.19.2 matplotlib:3.3.4 scipy:1.5.2")
    print("tqdm:4.61.2 ")
    print("===============================================\n")
    main()
    print("complete! [Enter > ]")
    input()