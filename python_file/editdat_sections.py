# python3.6.5 2022/10/11 yamashita

# SeamFEMの入力ファイル（sections作成）
import tkinter as tk
import tkinter.filedialog as fd
import numpy as np
import math
import os
from matplotlib import pyplot as plt


def input_dat():
    
    root = tk.Tk()
    root.withdraw()

    print("編集するda_ファイルを選択>>>")
    inf = fd.askopenfilename(
        filetypes=[("dat", ".da_")]
    )
    if inf == "":
        exit("[ERROR] No file selected")
    print("da_file: %s\n" %(inf))
    return inf

def input_vtx():
    root = tk.Tk()
    root.withdraw()

    print("入力する頂点データの選択>>>")
    inf = fd.askopenfilename(
        filetypes=[("xyz", "twodimension.xyz")]
    )
    if inf == "":
        exit("[ERROR] No file selected")
    print("vtx file: %s\n" %(inf))
    return inf

def output_dat():

    root = tk.Tk()
    root.withdraw()

    print("出力ファイルの作成")
    ouf = fd.asksaveasfilename(
        title = "名前を付けて保存",
        filetypes=[("dat", ".da_")],
        # initialdir=dpath,
        defaultextension = "da_"
    )
    if ouf == "":
        exit("[ERROR] No file selected")
    print("output da_file: %s\n" %(ouf))
    return ouf

def main():

    global MEAN, N

    MEAN = False #True -> 平均をとる　False -> 中央値をとる
    N = 50 #最大幅の分割数の指定

    inf = input_vtx()
    in_dat = input_dat()
    # ouf = output_dat()

    # 二次元頂点座標データの読みこみ
    with open(inf, "r") as f: 
        row = f.readlines()

    data = np.zeros(((int(len(row)/12), 12, 2))) #各断面ごとの頂点座標を格納
    for i in range(int(len(row)/12)):
        for j in range(12):
            rowdata = row[12*i + j].replace("\n", "")
            rowdata = rowdata.split(" ")
            data[i][j] = np.array([float(rowdata[0]), float(rowdata[1])])


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

    if MEAN == True:

        # 平均をとる
        uf = np.array([np.mean(uflange.T[0]), np.mean(uflange.T[1]), np.mean(uflange.T[2]), np.mean(uflange.T[3]), np.mean(uflange.T[4]), np.mean(uflange.T[5])])
        df = np.array([np.mean(dflange.T[0]), np.mean(dflange.T[1]), np.mean(dflange.T[2]), np.mean(dflange.T[3]), np.mean(dflange.T[4]), np.mean(dflange.T[5])])
        w = np.array([np.mean(web.T[0]), np.mean(web.T[1]), np.mean(web.T[2]), np.mean(web.T[3]), np.mean(web.T[4]), np.mean(web.T[5])])

    else:

        # 中央値をとる
        uf = np.array([np.median(uflange.T[0]), np.median(uflange.T[1]), np.median(uflange.T[2]), np.median(uflange.T[3]), np.median(uflange.T[4]), np.median(uflange.T[5])])
        df = np.array([np.median(dflange.T[0]), np.median(dflange.T[1]), np.median(dflange.T[2]), np.median(dflange.T[3]), np.median(dflange.T[4]), np.median(dflange.T[5])])
        w = np.array([np.median(web.T[0]), np.median(web.T[1]), np.median(web.T[2]), np.median(web.T[3]), np.median(web.T[4]), np.median(web.T[5])])

        
    # 分割した際に正方形状になるために縦横比から分割数を決定する。
    argm_data = np.array([uf, df, w])

    
    nn = np.zeros((3, 2)) #分割数を格納するための配列(0: 厚さ、1:幅)
    for i in range(3):
        var = np.mean(argm_data[i][-2:])
        hor = np.linalg.norm(argm_data[i][0:2] - argm_data[i][2:4])
        if var < hor:
            rat = var / hor
            nn[i] = np.array([math.ceil(N*rat), N])
        else:
            rat = hor / var
            nn[i] = np.array([N, math.ceil(N*rat)])
    

    # 編集datファイルのインプット
    with open(in_dat, "r") as f:
        row_dat = f.read()
        
    row_dat = row_dat.split("/sections")
    row_dat[1] = row_dat[1].split("/quake")[1]


    # 書き込み
    ouf_path = os.path.splitext(in_dat)[0] + "_editedSections.da_"
    with open(ouf_path, "w") as f:
        f.write(row_dat[0])
        f.write("/sections\n")
        print("{:10d}{:10d}{:10d}{:10d}{:15}{:10}{:10d}{:10d}{:10d}{:10d}{:10d}{:10d}{:10d}{:10d}{:10d}".format(1, 1, 3, 0, "", "", 0, 0, 1, 0, 0, 0, 0, 0, 0), file = f)
    
        for i in range(3):
            print("{:10d}{:10d}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10d}{:10d}" .format(int(nn[i][0]), int(nn[i][1]), argm_data[i][0], argm_data[i][1], argm_data[i][2], argm_data[i][3], argm_data[i][4], argm_data[i][5], 0, 0) ,file = f)
        
        f.write("/quake\n")
        f.write(row_dat[1])  
    
    print("complete! [Enter >]", end="")
    input()
    exit()


main()