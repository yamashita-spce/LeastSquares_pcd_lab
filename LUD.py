#最小二乗法による平面の検出に必要な数学的な処理
#2021/07

import numpy as np

#LU分解（A=LUとなるような、左下三角行列Lと右上三角行列Uを求める）
#Aは正方行列
def lud(A):

    #条件
    if len(A) != len(A[0]):
        print("[ERROR] LUD.py\nThe argument must be a square matrix")
        exit()

    #配列と変数の初期化
    ROW = len(A)
    L = np.zeros((ROW,ROW))
    U = np.zeros((ROW,ROW))
    Ai = np.zeros((ROW,ROW))
    Ai = np.copy(A)

    for i in range(ROW):

        #LUに値を代入
        for k in range(ROW - i):
            L[ROW - (k + 1)][i] = Ai[ROW - (k + 1)][i]
            U[i][ROW - (k + 1)] = Ai[i][ROW - (k + 1)]/Ai[i][i]
        
            if i == (ROW - (k + 1)):
                U[i][ROW - (k + 1)] = 1
        
        #終了条件
        if (ROW - (i + 1)) == 0:
            break
        
        #Ai <-- Ai - l(1)*u(0)_T の変換はAi配列に直接変換値を代入している
        for j in range(ROW - (i + 1)):
            for m in range(ROW - (i + 1)):
                #Ai <-- Ai - l_1*u0_T
                Ai[ROW - (1 + j)][ROW - (1 + m)] -=  L[ROW - (j + 1)][i]* U[i][ROW - (m + 1)]

    #戻り値はLとU（numpy二次元配列）
    return(L, U)


# LUX = Yを満たすXの解を求める。
# ただしL,UはLU分解されたn次元numpy正方行列のみ引数にとる。Yは1次元numpy行列をとる
def lusolv(L, U, Y):
    
    #配列、変数の初期化
    ROW = len(L)
    B = np.zeros(ROW)
    X = np.zeros(ROW)

    #LB = YからBを求める
    for i in range(ROW):
    
        B[i] = Y[i]/L[i][i]
        if i == 0:
            continue
        
        for j in range(i):
            B[i] -= (L[i][i - j - 1]/L[i][i])*B[i - j - 1]

    # UX = BからXを求める
    X = np.zeros(ROW)

    for i in range(ROW):
        j = ROW - 1 - i
    
        X[j] = B[j]
        if i == 0:
            continue
    
        for k in range(i):
            X[j] -= U[j][ROW - (k + 1)]*X[ROW - (k + 1)]

    #戻り値は解となる一次元numpy行列
    return(X)

# 3次元点群から統計データを返す関数。
# 引数xyzはnumpy二次元配列
def datst_xyz(xyz):

    #変数、配列の初期化
    n = len(xyz)
    sxyz = np.zeros(3)
    sxyz2 = np.zeros(3)
    sxyyz = np.zeros(3)

    for xyz_i in xyz:

        sxyz[0] += xyz_i[0]
        sxyz[1] += xyz_i[1]
        sxyz[2] += xyz_i[2]
        sxyz2[0] += xyz_i[0]*xyz_i[0]
        sxyz2[1] += xyz_i[1]*xyz_i[1]
        sxyz2[2] += xyz_i[2]*xyz_i[2]
        sxyyz[0] += xyz_i[0]*xyz_i[1]
        sxyyz[1] += xyz_i[1]*xyz_i[2]
        sxyyz[2] += xyz_i[2]*xyz_i[0]

    # nは総データ数、
    # sxyzは各座標の合計値、sxyz2は各座標の2乗の合計値、
    # sxyyzは各座標のxy、yz、zxの積の合計値
    return(n, sxyz, sxyz2, sxyyz)


# 任意の軸に対して最小二乗平面による連立方程式の係数正方行列Aと縦行列Yを求める
# 引数はLUD.datst_xyzの戻り値（n, sxyz, sxyz2, sxyyz）に対応.SHAFTは軸(x,y,z)の指定
def lstsqrpln(n, sxyz, sxyz2, sxyyz, SHAFT):

    #軸分岐
    if SHAFT == "x":
        #a+by+cz = x --> AX = Y
        A = np.array([[n,       sxyz[1],  sxyz[2]],
                      [sxyz[1], sxyz2[1], sxyyz[1]],
                      [sxyz[2], sxyyz[1], sxyz2[2]]])

        Y = np.array([sxyz[0], sxyyz[0], sxyyz[2]])
    
    elif SHAFT == "y":
        #a+bx+cz = y --> AX = Y
        A =  np.array([[n,       sxyz[0],  sxyz[2]],
                      [sxyz[0], sxyz2[0], sxyyz[2]],
                      [sxyz[2], sxyyz[2], sxyz2[2]]])
        
        Y = np.array([sxyz[1], sxyyz[0], sxyyz[1]])

    elif SHAFT == "z":
        #a+bx+cy = z --> AX = Y
        A = np.array([[n,       sxyz[0],  sxyz[1]],
                      [sxyz[0], sxyz2[0], sxyyz[0]],
                      [sxyz[1], sxyyz[0], sxyz2[1]]])
        
        Y = np.array([sxyz[2], sxyyz[2], sxyyz[1]])

    else:
        print("[ERORR] LUD.lstsqrpln()\nArgument SHAFT is an unknown value")
        exit()

    #戻り値Aは3行正方numpy二次元配列、Yは3行numpy一次元配列
    return(A, Y)
