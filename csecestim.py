#断面検出
#python 環境一覧（2021/07)
# python : version = 3.6.5  Anaconda : version = 4.10.3
# Library : numpy=1.19.2, tkinker=8.6.10, matplotlib=3.3.4, 
# local module : LUD (LU分解、最小二乗平面）, pathio (tk)

import numpy as np
from matplotlib import markers, pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#file ioモジュール
import pathio
#LU分解モジュール
import LUD 

#点群ファイルを選択し、初期表示
def io_begin_plot():

    #データをインポート
    inf = pathio.infle()
    print("\nSelected : %s\n" %(inf))
    with open(inf, "r") as f:
        dat = f.readlines()

    le = len(dat)
    xyz = np.zeros((le,3))

    for i,li in enumerate(dat):
        li = li.replace(" \n","")
        data = li.split(" ")
        a,b,c = float(data[0]),  float(data[1]),  float(data[2])
        xyz[i] = [a,b,c]

    # フォントの種類とサイズを設定
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'Times New Roman'

    # グラフの入れ物を用意
    fig = plt.figure()
    ax = Axes3D(fig)

    # 軸のラベルを設定
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    

    # データプロット
    xyz_t = xyz.T
    ax.scatter3D(xyz_t[0], xyz_t[1], xyz_t[2], marker=".",label='Dataset')
    plt.legend()
    plt.show()
    plt.close()

    #検出したい平面の入力
    while True:
        print("Reference plane input [xy:yz:zx] >")
        rp = input()
        if (rp != "xy" and rp != "yz" and rp != "zx"):
            print("[ERROR] Enter one of xy, yz, xz\n")
            continue
        elif rp == "xy":
            rs = "z"
            break
        elif rp == "yz":
            rs = "x"
            break
        elif rp == "zx":
            rs = "y"
            break

    return(rp, rs, xyz, xyz_t)

#断面を構成する点群を検出
def crosecdetct(rp, xyz, xyz_t):

    #基準平面の指定
    PLN, PLNS = 0, ""
    if rp == "xy":
        PLN, PLNS = 2, "z"
    elif rp == "yz":
        PLN, PLNS = 0, "x"
    elif rp == "zx":
        PLN, PLNS = 1, "y"
    
    #精度の指定（デフォルト2000) <-- 点群密度によって値を変える必要あり
    d = 2000

    #点群がとる座標の最大値及び最小値
    pmin, pmax = min(xyz_t[PLN]), max(xyz_t[PLN])

    #座標における点密度の平均値
    page = int(len(xyz)/d)

    print("\nDetected point cloud data >>")
    print("\nMaximum value of %s : %f" %(PLNS, pmax))
    print("Minimum value of %s : %f" %(PLNS, pmin))
    print("Average density of %s : %d [1/%d]" %(PLNS, page, d))


    #反復で用いる変数、配列の定義
    idx = pmin
    p_list0 = [] 
    pcd, ix = np.zeros(d + 2), np.zeros(d + 2)
    rge = float((pmax - pmin)/d)

    #座標が小さい面を構成している点群データのインデックスを取得
    k, kk = 0, 0
    for i in range(d + 2):
        
        for j,p in enumerate(xyz_t[PLN]):
            if ((idx - rge) <= p and p < idx):
                p_list0.append(j)
                pcd[i] += 1
            
        #点密度が平均を上回ったときにkk=1とする
        if pcd[i] > page:
            kk = 1
        #点密度が平均を下回った回数をkにキャッシュ
        elif (pcd[i] < page and kk == 1):
            k += 1
            #k > d*(1/50)となったら終了(最適化)　<-- 点群の状態により変更する必要あり
            if k > d*(1/50):
                break
        
        ix[i] = idx
        idx += rge

    if DETAIL == True:

        #検出した点群密度をグラフ表示
        fig0, ax0 = plt.subplots()
        ax0.plot(ix, pcd)
        ax0.set_xlim(pmin, pmin + d*rge/20)
        ax0.set_xlabel(PLNS)
        ax0.set_ylabel('Density distribution')
        
        ax0.grid()
        plt.show()
        plt.close()

    #検出した点群インデックス（p_list0)をつかって、新たな点群配列を作成
    points0 = np.zeros((len(p_list0),3))

    for i,idx in enumerate(p_list0):
        points0[i] = xyz[idx]

    if DETAIL == True:

        # グラフの入れ物を用意する。
        fig1 = plt.figure()
        ax1 = Axes3D(fig1)

        # 軸のラベルを設定する。
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_zlabel('z')

        # データプロットする。
        points0_t = points0.T
        ax1.scatter3D(points0_t[0], points0_t[1], points0_t[2] ,marker=".", label='Dataset')
        plt.legend()
        plt.show()
        plt.close()


    #反復で用いる変数、配列の定義、初期化
    p_list1 = []
    idx = pmax
    pcd, ix = np.zeros(d + 2), np.zeros(d + 2)

    #  ＜座標が大きい面を構成している点群データのインデックスを取得＞
    k, kk = 0, 0
    for i in range(d + 2):
        
        for j,p in enumerate(xyz_t[PLN]):
            if (idx < p and p <= idx + rge):
                p_list1.append(j)
                pcd[d + 1 - i] += 1
            
        #点密度が平均を上回ったときにkk=1とする
        if pcd[d + 1 - i] > page:
            kk = 1
        #点密度が平均を下回った回数をkにキャッシュ
        elif (pcd[d + 1 - i] < page and kk == 1):
            k += 1
            #k > d*(1/50)となったら終了(最適化)
            if k > d*(1/50):
                break
        
        ix[d + 1 - i] = idx
        idx -= rge

    if DETAIL == True:

        #検出した点群密度をグラフ表示
        fig2, ax2 = plt.subplots()
        ax2.plot(ix, pcd)
        ax2.set_xlim(pmax - d*rge/20, pmax)
        ax2.set_xlabel(PLNS)
        ax2.set_ylabel('Density distribution')
    
        ax2.grid()
        plt.show()
        plt.close()

    #検出した点群インデックス（p_list0)をつかって、新たな点群配列を作成
    points1 = np.zeros((len(p_list1),3))

    for i,idx in enumerate(p_list1):
        points1[i] = xyz[idx]

    if DETAIL == True:

        fig3 = plt.figure()
        ax3 = Axes3D(fig3)
        ax3.set_xlabel('x')
        ax3.set_ylabel('y')
        ax3.set_zlabel('z')

        points1_t = points1.T
        ax3.scatter3D(points1_t[0], points1_t[1], points1_t[2] ,marker=".", label='Dataset')
        plt.legend()
        plt.show()
        plt.close()

    return(points0, points1)

#平面に点をそらせる
def warptosurface(xyzw, X, rsw):
    
    xyzk = np.zeros((len(xyzw), 3))

    #基準軸で分岐
    if rsw == "x":
        for i,p in enumerate(xyzw):
            xyzk[i][0] = X[0] + X[1]*p[1] + X[2]*p[2]
            xyzk[i][1] = p[1]
            xyzk[i][2] = p[2]

    elif rsw == "y":
        for i,p in enumerate(xyzw):
            xyzk[i][0] = p[0]
            xyzk[i][1] = X[0] + X[1]*p[0] + X[2]*p[2]
            xyzk[i][2] = p[2]
 
    elif rsw == "z":
        for i,p in enumerate(xyzw):
            xyzk[i][0] = p[0]
            xyzk[i][1] = p[1]
            xyzk[i][2] = X[0] + X[1]*p[0] + X[2]*p[1]
    
    return(xyzk)
    

# 断面を構成している点群の正規化
def normalofsecpcd(xyzn, X, rs):

    #正規化するために、点群を回転移動させる(複素平面を用いる)
    xyz_copy = np.copy(xyzn)
    new_xyz = np.zeros((len(xyzn),3))

    rad0 = np.pi/2 - np.abs(np.arctan(1/X[1]))
    rad1 = np.pi/2 - np.abs(np.arctan(1/X[2]))
    sina, cosa = np.sin(rad0), np.cos(rad0)
    sinb, cosb = np.sin(rad1), np.cos(rad1)

    #基準軸で分岐
    if rs == "x":

        for i,ex in enumerate(xyz_copy):
            #x,y,z = ex[0] - X[0], ex[1], ex[2] <-- 座標を(0, 0, 0)に正規化する
            x,y,z = ex[0], ex[1], ex[2]

            new_xyz[i][0] = (x*cosa - y*sina)*cosb - z*sinb
            new_xyz[i][1] = x*sina + y*cosa
            new_xyz[i][2] = (x*cosa - y*sina)*sinb + z*cosb
    
    if rs == "y":

        for i,ex in enumerate(xyz_copy):
            #x,y,z = ex[0], ex[1] - X[0], ex[2]
            x,y,z = ex[0], ex[1], ex[2]

            new_xyz[i][0] = y*sina + x*cosa
            new_xyz[i][1] = (y*cosa - x*sina)*cosb - z*sinb
            new_xyz[i][2] = (y*cosa - x*sina)*sinb + z*cosb

    if rs == "z":

        for i,ex in enumerate(xyz_copy):
            #x,y,z = ex[0], ex[1], ex[2] - X[0]
            x,y,z = ex[0], ex[1], ex[2]

            new_xyz[i][0] = z*sina + x*cosa
            new_xyz[i][1] = (z*cosa - x*sina)*sinb + y*cosb
            new_xyz[i][2] = (z*cosa - x*sina)*cosb - y*sinb

    return(new_xyz, X[0])

#コメントを返す関数
def comment01(X):
    
    print("[WARNING] Normalization creates a large error")
    print(">> Still do normalization [y:n] >")
    while True:
        anser01 = input()
        if anser01 == "y":
            X[0],X[1],X[2] = 0, 0, 0
            break
        if anser01 == "n":
            break
        else:
            print("[ERROR] Enter one of y, n\n")
            continue
    
    return(X)

#平面を推定し、正規化した同一平面上の点群を返す関数
def subesurface(dpcd, dsm):

    # 配列の作成
    n0, sxyz0, sxyz20, sxyyz0 = LUD.datst_xyz(dpcd)
    A, Y = LUD.lstsqrpln(n0, sxyz0, sxyz20, sxyyz0, dsm)
    
    # LU分解から解を求める
    L, U = LUD.lud(A)
    X = LUD.lusolv(L, U, Y)
    #warptosurface(dpcd, X, dsm)

    #検出した平面の方程式を出力
    print("\nDetection section equation >>\n ")
    if dsm == "x":
        print("x = %f + %fy + %fz" %(X[0], X[1], X[2]))
    elif dsm == "y":
        print("y = %f + %fx + %fz" %(X[0], X[1], X[2]))
    elif dsm == "z":
        print("z = %f + %fx + %fy" %(X[0], X[1], X[2]))

    # 点群の正規化
    dpcd_new, adder = normalofsecpcd(dpcd, X, dsm)

    # 正規化した点群について平面を検出
    n0, sxyz0, sxyz20, sxyyz0 = LUD.datst_xyz(dpcd_new)
    A, Y = LUD.lstsqrpln(n0, sxyz0, sxyz20, sxyyz0, dsm)
    L, U = LUD.lud(A)
    X = LUD.lusolv(L, U, Y)
    
    print("                  _")
    print("                 | |")
    print("Normalization    | |")
    print("                \   /")
    print("                 \ /\n")

    if dsm == "x":
        print("x = %f + %fy + %fz" %(X[0], X[1], X[2]))
    elif dsm == "y":
        print("y = %f + %fx + %fz" %(X[0], X[1], X[2]))
    elif dsm == "z":
        print("z = %f + %fx + %fy" %(X[0], X[1], X[2]))


    # 係数を0に正規化する

    #断面の座標を(0, 0, 0)に正規化
    #if (X[0] <= 0.001) and (X[1] <= 0.001) and (X[2] <= 0.001):
    #   X[0], X[1], X[2] = 0, 0, 0

    #断面の基準座標は変更しない
    if (X[1] <= 0.001) and (X[2] <= 0.001):
        X[1], X[2] = 0, 0
        print("\nSuccessful normalization")

    else:
        X = comment01(X)
        
    finpcd = warptosurface(dpcd_new, X, dsm)

    return(finpcd, adder)

#3次元の点群をグラフで表示する
def matplot(p):
    
    points_t = p.T

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    ax.scatter3D(points_t[0], points_t[1], points_t[2] ,marker=".", label='Dataset')
    plt.legend()
    plt.show()
    plt.close()

def main():
    print("========================================================================")
    print("CSECESTIM.EXE python : version = 3.6.5  Anaconda : version = 4.10.3")
    print("Library : numpy=1.19.2, tkinker=8.6.10, matplotlib, mpl_toolkits")
    print("local module : LUD (LU分解、最小二乗平面）, pathio (tk)")
    print("========================================================================")
    
    global DETAIL 

    #DETAIL = False
    DETAIL = True

    #座標データと検出したい断面を返す関数
    rpm, rsm, xyzm, xyz_tm = io_begin_plot()

    # 断面を構成する点群を切り取る
    cs0, cs1 = crosecdetct(rpm, xyzm, xyz_tm)

    # 点群が構成している平面を推定する

    sfpcd0, ditnce0 = subesurface(cs0, rsm)
    matplot(sfpcd0)
    sfpcd1, ditnce1 = subesurface(cs1, rsm)
    matplot(sfpcd1)

if __name__ == "__main__":
    main()

