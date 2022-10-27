# SeamFEMの入力ファイル（新規作成, 編集）

import tkinter as tk
import tkinter.filedialog as fd
import numpy as np
from matplotlib import pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D

# 特徴点座標ファイルの入力
def input_vtx():
    global INF, DIR

    print("choose input tmp file ", end=">  ")
    root = tk.Tk()
    root.withdraw()
    INF = fd.askopenfilename(
        filetypes=[("tmp", "twodimension.tmp")]
    )
    DIR = os.path.split(INF)[0] + "/"
    print("%s" %(INF))

    try: 
        with open(INF, "r") as f:
            row = f.readlines()
    except:
        exit("[ERROR] file %s doesn't exist." %(INF))
    
    data = np.zeros(((int(len(row)/13), 12, 2))) #各断面ごとの頂点座標を格納
    G_data = np.zeros((int(len(row)/13), 3)) #三次元重心データの格納

    for i in range(int(len(row)/13)):
        rowdata = row[13*i].replace("\n", "")
        rowdata = rowdata.split(" ")
        G_data[i] = np.array([float(rowdata[0]), float(rowdata[1]), float(rowdata[2])])
        
        for j in range(12):
            rowdata = row[13*i + j + 1].replace("\n", "")
            rowdata = rowdata.split(" ")
            data[i][j] = np.array([float(rowdata[0]), float(rowdata[1])])

    return data, G_data

    
# datファイルの選択（編集するなら）
def input_dat():
    global EDIT

    print("choose input dat file ", end=">  ")
    root = tk.Tk()
    root.withdraw()
    EDIT = fd.askopenfilename(
        filetypes=[("dat", ".da_")]
    )
    if EDIT != "":
        print("%s" %(EDIT))
    else:
        print("新規作成")

def select_conf():
    global CS

    while(1):
        print("\n断面構成を指定 [0:平均値, 1:中央値, 2:断面を統一しない]", end=" >  ")
        CS = input()
        if CS=="0" or CS=="1" or CS=="2":
            break
        else:
            print("[WARNING] 0, 1, 2 のいずれかを入力")
    CS = int(CS)

def main():

    global EDIT, INF, DIR, N, CS, DESCRIPSION


    # ======================================
    N = 50 #最大幅の分割数の指定
    DESCRIPSION = False #節点の描写
    # ====================================== 


    data, G_data = input_vtx()
    input_dat()
    select_conf()

    # 座標系を各重心点（0,0）に合わせる
    for i in range(len(data)):
        mean = np.array([np.mean(data[i][8:].T[0]), np.mean(data[i][8:].T[1])])
        
        for j in range(12):
            data[i][j] -= mean


    # 各面を構成する長方形の厚さとその端中心座標を求める
    uflange = np.zeros((len(data), 6)) #カラム補足　ー＞　0,1: 左板中心座標(x,y) 、2.3:右板中心座標(x.,y)、4:左板厚、5:右板厚
    dflange = np.zeros((len(data), 6))
    web = np.zeros((len(data), 6))

    for i in range(len(data)):
        ufl_mean = np.array([np.mean(data[i][0:2].T[0]), np.mean(data[i][0:2].T[1])])
        ufr_mean = np.array([np.mean(data[i][2:4].T[0]), np.mean(data[i][2:4].T[1])])
        dfl_mean = np.array([np.mean(data[i][4:6].T[0]), np.mean(data[i][4:6].T[1])])
        dfr_mean = np.array([np.mean(data[i][6:8].T[0]), np.mean(data[i][6:8].T[1])])
        lw_mean = np.array([np.mean(data[i][[8,10]].T[0]), np.mean(data[i][[8,10]].T[1])])
        rw_mean = np.array([np.mean(data[i][[9,11]].T[0]), np.mean(data[i][[9,11]].T[1])])
        ufl_lig = np.linalg.norm(data[i][0] - data[i][1])
        ufr_lig = np.linalg.norm(data[i][2] - data[i][3])
        dfl_lig = np.linalg.norm(data[i][4] - data[i][5])
        dfr_lig = np.linalg.norm(data[i][6] - data[i][7])
        lw_lig = np.linalg.norm(data[i][8] - data[i][10])
        rw_lig = np.linalg.norm(data[i][9] - data[i][11])
        
        uflange[i] = np.array([ufl_mean[0], ufl_mean[1], ufr_mean[0], ufr_mean[1], ufl_lig, ufr_lig])
        dflange[i] = np.array([dfl_mean[0], dfl_mean[1], dfr_mean[0], dfr_mean[1], dfl_lig, dfr_lig])
        web[i] = np.array([lw_mean[0], lw_mean[1], rw_mean[0], rw_mean[1], lw_lig, rw_lig])

    if CS == 0:
        # 平均をとる
        uf = np.array([np.mean(uflange.T[0]), np.mean(uflange.T[1]), np.mean(uflange.T[2]), np.mean(uflange.T[3]), np.mean(uflange.T[4]), np.mean(uflange.T[5])])
        df = np.array([np.mean(dflange.T[0]), np.mean(dflange.T[1]), np.mean(dflange.T[2]), np.mean(dflange.T[3]), np.mean(dflange.T[4]), np.mean(dflange.T[5])])
        w = np.array([np.mean(web.T[0]), np.mean(web.T[1]), np.mean(web.T[2]), np.mean(web.T[3]), np.mean(web.T[4]), np.mean(web.T[5])])
    elif CS == 1:
        # 中央値をとる
        uf = np.array([np.median(uflange.T[0]), np.median(uflange.T[1]), np.median(uflange.T[2]), np.median(uflange.T[3]), np.median(uflange.T[4]), np.median(uflange.T[5])])
        df = np.array([np.median(dflange.T[0]), np.median(dflange.T[1]), np.median(dflange.T[2]), np.median(dflange.T[3]), np.median(dflange.T[4]), np.median(dflange.T[5])])
        w = np.array([np.median(web.T[0]), np.median(web.T[1]), np.median(web.T[2]), np.median(web.T[3]), np.median(web.T[4]), np.median(web.T[5])])


    # 分割した際に正方形状になるために縦横比から分割数を決定する。
    if CS==0 or CS==1:
        argm_data = np.zeros((3, 6))
        argm_data = np.array([uf, df, w])
        nn = np.zeros((3, 2)) #分割数を格納するための配列(0: 厚さ、1:幅)

        for i in range(3):
            var = np.mean(argm_data[i][-2:])
            hor = np.linalg.norm(argm_data[i][0:2] - argm_data[i][2:4])
            if var < hor:
                rat = var / hor
                nn[i] = np.array([int(N*rat), N])
            else:
                rat = hor / var
                nn[i] = np.array([N, int(N*rat)])
    
    elif CS==2:
        argms_data = np.zeros((len(data), 3, 6))
        for i in range(len(data)):
            argms_data[i] = np.array([uflange[i], dflange[i], web[i]])

        nnn = np.zeros(((len(data), 3, 2))) #分割数を格納するための配列(0: 各断面, 1: 厚さ、2:幅)

        for i in range(len(data)):
            for j in range(3):
                var = np.mean(argms_data[i][j][-2:])
                hor = np.linalg.norm(argms_data[i][j][0:2] - argms_data[i][j][2:4])
                if var < hor:
                    rat = var / hor
                    nnn[i][j] = np.array([int(N*rat), N])
                else:
                    rat = hor / var
                    nnn[i][j] = np.array([N, int(N*rat)])
            

    # 節点作成
    node_xyz = np.zeros((len(G_data) + 1, 3))

    # 各隣同士の点の中心地と両端の単位ベクトルを求める
    center = np.zeros((len(G_data) - 1, 3)) #x, y, z, 
    vector = np.zeros((2, 3)) #[x, y, z]単位ベクトル

    for i in range(len(G_data) - 1):
        center[i] = np.array([np.mean(G_data.T[0][i:i+2]), np.mean(G_data.T[1][i:i+2]), np.mean(G_data.T[2][i:i+2])])
        
    vector[0] = (G_data[1] - G_data[0]) / 2
    vector[1] = (G_data[-1] - G_data[-2]) / 2

    node_xyz[0] = G_data[0] - vector[0]
    node_xyz[-1] = G_data[-1] + vector[1]
    node_xyz[1:-1] = center

    # 節点と重心点の描写
    if DESCRIPSION:
        print("plt 3Dplot...")
        fig = plt.figure()
        ax1 = Axes3D(fig)
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_zlabel('z')
        ax1.scatter3D(node_xyz.T[0], node_xyz.T[1], node_xyz.T[2], label="nodes")
        ax1.scatter3D(G_data.T[0], G_data.T[1], G_data.T[2], label="center of gravity")
        plt.legend()
        plt.show()


    # 新規書き込み
    if EDIT == "":
        if CS==0 or CS==1:
            print("新規ファイル: %s" %(DIR + "create.da_"))
            with open(DIR + "create.da_", "w") as f:
                print("/nodes", file=f)
                for i in range(len(node_xyz)):
                    print("{:10d}{:10d}{:10d}{:>15f}{:>15f}{:>15f}{:10d}{:10d}{:10d}".format(1+i, 0, 0, node_xyz[i][0], node_xyz[i][1], node_xyz[i][2], 0, 0, 0), file=f)
                
                print("/elements", file=f)
                for i in range(len(node_xyz) - 1):
                    print("{:10d}{:10d}{:10d}{:10d}{:10d}{:10d}".format(1+i, 0, 0, 1, 1+i, 2+i), file=f)
            
                print("/sections", file=f)
                print("{:>10d}{:>10d}{:>10d}{:>10d}{:>15}{:>10}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}".format(1, 1, 3, 0, "", "", 0, 0, 1, 0, 0, 0, 0, 0, 0), file = f)
            
                for i in range(3):
                    print("{:>10d}{:>10d}{:>10f}{:>10f}{:>10f}{:>10f}{:>10f}{:>10f}{:>10d}{:>10d}" .format(int(nn[i][0]), int(nn[i][1]), argm_data[i][0], argm_data[i][1], argm_data[i][2], argm_data[i][3], argm_data[i][4], argm_data[i][5], 0, 0) ,file = f)
        elif CS==2:
            print("新規ファイル: %s" %(DIR + "create_dif.da_"))
            with open(DIR + "create_dif.da_", "w") as f:
                print("/nodes", file=f)
                for i in range(len(node_xyz)):
                    print("{:10d}{:10d}{:10d}{:>15f}{:>15f}{:>15f}{:10d}{:10d}{:10d}".format(1+i, 0, 0, node_xyz[i][0], node_xyz[i][1], node_xyz[i][2], 0, 0, 0), file=f)
                
                print("/elements", file=f)
                for i in range(len(node_xyz) - 1):
                    print("{:10d}{:10d}{:10d}{:10d}{:10d}{:10d}".format(1+i, 0, 0, 1+i, 1+i, 2+i), file=f)
                
                for i in range(len(data)):
                    print("/sections", file=f)
                    print("{:>10d}{:>10d}{:>10d}{:>10d}{:>15}{:>10}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}".format(i+1, 1, 3, 0, "", "", 0, 0, 1, 0, 0, 0, 0, 0, 0), file = f)
                    for j in range(3):
                        print("{:>10d}{:>10d}{:>10f}{:>10f}{:>10f}{:>10f}{:>10f}{:>10f}{:>10d}{:>10d}" .format(int(nnn[i][j][0]), int(nnn[i][j][1]), argms_data[i][j][0], argms_data[i][j][1], argms_data[i][j][2], argms_data[i][j][3], argms_data[i][j][4], argms_data[i][j][5], 0, 0) ,file=f)
    
    # 上書き保存  
    else:
        with open(EDIT, "r") as f:
            row_dat = f.read()
        
        try:
            row_tmp, row_lmp = row_dat.split("/nodes")
            row_lmp = row_lmp.split("/")[1:]
            row_lmp = "/" + "/".join(row_lmp)
            row_dat = row_tmp + row_lmp
        except:
            print("[WARNING] /nodesを新規作成します.")
                
        try:
            row_tmp, row_lmp = row_dat.split("/elements")
            row_lmp = row_lmp.split("/")[1:]
            row_lmp = "/" + "/".join(row_lmp)
            row_dat = row_tmp + row_lmp
        except:
            print("[WARNING] /elementsを新規作成します.")

        try:
            row_tmp, row_lmp = row_dat.split("/sections")
            row_lmp = row_lmp.split("/")[1:]
            row_lmp = "/" + "/".join(row_lmp)
            row_dat = row_tmp + row_lmp
        except:
            print("[WARNING] /sectionsを新規作成します.")
        
        edir = os.path.splitext(EDIT)[0]

        if CS==0 or CS==1:
            print("上書きファイル: %s" %(edir + "_edit.da_"))
            with open(edir + "_edit.da_", "w") as f:
                print("/nodes", file=f)
                for i in range(len(node_xyz)):
                    print("{:10d}{:10d}{:10d}{:>15f}{:>15f}{:>15f}{:10d}{:10d}{:10d}".format(1+i, 0, 0, node_xyz[i][0], node_xyz[i][1], node_xyz[i][2], 0, 0, 0), file=f)
                
                print("/elements", file=f)
                for i in range(len(node_xyz) - 1):
                    print("{:10d}{:10d}{:10d}{:10d}{:10d}{:10d}".format(1+i, 0, 0, 1, 1+i, 2+i), file=f)

                print("/sections", file=f)
                print("{:>10d}{:>10d}{:>10d}{:>10d}{:>15}{:>10}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}".format(1, 1, 3, 0, "", "", 0, 0, 1, 0, 0, 0, 0, 0, 0), file = f)
            
                for i in range(3):
                    print("{:>10d}{:>10d}{:>10f}{:>10f}{:>10f}{:>10f}{:>10f}{:>10f}{:>10d}{:>10d}" .format(int(nn[i][0]), int(nn[i][1]), argm_data[i][0], argm_data[i][1], argm_data[i][2], argm_data[i][3], argm_data[i][4], argm_data[i][5], 0, 0) ,file = f)
        
                print(row_dat, file=f)

        elif CS==2:
            print("上書きファイル: %s" %(edir + "_edit_dif.da_"))
            with open(edir + "_edit_dif.da_", "w") as f:
                print("/nodes", file=f)
                for i in range(len(node_xyz)):
                    print("{:10d}{:10d}{:10d}{:>15f}{:>15f}{:>15f}{:10d}{:10d}{:10d}".format(1+i, 0, 0, node_xyz[i][0], node_xyz[i][1], node_xyz[i][2], 0, 0, 0), file=f)
                
                print("/elements", file=f)
                for i in range(len(node_xyz) - 1):
                    print("{:10d}{:10d}{:10d}{:10d}{:10d}{:10d}".format(1+i, 0, 0, 1+i, 1+i, 2+i), file=f)
                
                for i in range(len(data)):
                    print("/sections", file=f)
                    print("{:>10d}{:>10d}{:>10d}{:>10d}{:>15}{:>10}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}{:>10d}".format(i+1, 1, 3, 0, "", "", 0, 0, 1, 0, 0, 0, 0, 0, 0), file = f)
                    for j in range(3):
                        print("{:>10d}{:>10d}{:>10f}{:>10f}{:>10f}{:>10f}{:>10f}{:>10f}{:>10d}{:>10d}" .format(int(nnn[i][j][0]), int(nnn[i][j][1]), argms_data[i][j][0], argms_data[i][j][1], argms_data[i][j][2], argms_data[i][j][3], argms_data[i][j][4], argms_data[i][j][5], 0, 0) ,file=f)
                
                print(row_dat, file=f)

main()
print("\ncomplete! ", end="[Enter >]")
input()

    

    

