#完全点群　extension

#詳細：
from scipy.optimize import minimize
import numpy as np
import math
import matplotlib.pyplot as plt


#断面点群の形状判別を行うクラス
class typeInference:

    #断面二次元点群
    pcd = np.zeros((0, 2))

    def __init__(self, p):
        self.pcd = p
    
    # Ｈ形のフィッティング->Ｔ型のフィッティングを行い,目的関数の返り値（fittingの良さを表すパラメータ値)が110％以上大きくなった場合、H型と推定する
    #T型なら0 H型なら1を返す
    def discrimination(self):

        hParameterCount = 6
        tParameterCount = 4

        hidx = self.HEvaluationValue()
        tidx = self.TEvaluationValue()
        n = len(self.pcd) #点数

        # パラメータ考慮した誤差値
        hfn = math.sqrt(hidx*hParameterCount / n)
        tfn = math.sqrt(tidx*tParameterCount / n)

        print("\n[返り値] H -> %f,  T -> %f" %(hidx, tidx))
        print("点群総数: %d 点, 標準偏差: H -> %f, T -> %f" %(n, math.sqrt(hidx / n), math.sqrt(tidx / n)))
        print("parameter数を考慮した誤差値: H -> %f, T -> %f" %(hfn, tfn))

        if hfn < tfn:
            return 1
        else:
            return 0


    #H型の評価を返す関数
    def HEvaluationValue(self):

        #H型初期値の設定
        #p = [Gx, Gy, h, rad] pxy = [x, y]
        p = np.array([np.average(self.pcd.T[0]), np.average(self.pcd.T[1]), np.abs(np.max(self.pcd.T[1]) - np.min(self.pcd.T[1])), 0, 0, 0])
        bo = np.abs(np.max(self.pcd.T[0]) - np.min(self.pcd.T[0]))

        P = minimize(self.Hobj_func, p, args=[self.pcd, bo], method="powell")
        print("H:")
        print(P)
     
        return P.fun


    #T型の評価を返す関数
    def TEvaluationValue(self):

        #T型初期値の設定
        #p = [Cx, Cy, Θ]
        p = np.array([np.mean(self.pcd.T[0]), np.max(self.pcd.T[1]), 0, 0])
        P = minimize(self.Tobj_func, p, args=self.pcd, method="powell")
        print("T:")
        print(P)

        return P.fun


    def Hobj_func(self, p, pxy): #H型目的関数
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
    
    #T型目的関数
    def Tobj_func(self, p, xy): 
    #パラメータ p = [Cx, Cy, α, β] xy = [x, y]
        SS = 0
        
        #直線の係数(ax + by + c = 0)
        tana = np.tan(p[2])
        tanb = np.tan(p[3])
        I = np.array([-1, tana, p[0] - tana*p[1]])
        II = np.array([tanb, -1, p[1] - tanb*p[0]])
        
        
        for i in range(len(xy)):
            
            ss = np.zeros(2)
            #点と直線の距離を求める
            ss[0] = (I[0]*xy[i][0] + I[1]*xy[i][1] + I[2])**2 / (I[0]**2 + I[1]**2)
            ss[1] = (II[0]*xy[i][0] + II[1]*xy[i][1] + II[2])**2 / (II[0]**2 + II[1]**2)
            
            SS += np.min(ss)
            
        return SS





# T型形状の処理オブジェクト
class TShapeProcessing:

    #二次元点群
    p = np.zeros((0, 2))
    #点情報(x, yの各々の最大最小値を格納)
    cfg = np.zeros((2, 2))
    #端子点
    A, B, C, D = np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2)

    #推定線分の係数行列
    I = np.zeros((2, 3))

    #点群グループ
    Ip = np.zeros((0, 2)) #web要素となる点群
    IIp = np.zeros((0, 2)) #フランジ要素となる点群

    Ap = np.zeros((0, 2)) #フランジA（上側）(四分位範囲)
    Bp = np.zeros((0, 2)) #フランジB（下側）(四分位範囲)
    Cp = np.zeros((0, 2)) #webC（左側）(四分位範囲)
    Dp = np.zeros((0, 2)) #webD（右側）(四分位範囲)

    # 特徴点
    Au = np.zeros(2) #フランジ左上
    Ad = np.zeros(2) #フランジ左下
    Bu = np.zeros(2) #フランジ右上
    Bd = np.zeros(2) #フランジ右下
    Cl = np.zeros(2) #フランジ-ウェブ接合点左
    Cr = np.zeros(2) #フランジ-ウェブ接合点右
    Dl = np.zeros(2) #ウェブ左下
    Dr = np.zeros(2) #ウェブ右下


    def __init__(self, p):
        self.p = p
        self.cfg = np.array([[np.max(p.T[0]), np.min(p.T[0])],[np.max(p.T[1]), np.min(p.T[1])]])

    #webの中心座標を返す
    def webCP(self, DISPLAY=False):

        xy = np.zeros(2)
        xy[0] = np.mean(np.array([self.Cl[0], self.Cr[0], self.Dl[0], self.Dr[0]]))
        xy[1] = np.mean(np.array([self.Cl[1], self.Cr[1], self.Dl[1], self.Dr[1]]))

        return xy

    #輪郭点群を返す
    def contourPcd(self, N=100, DISPLAY=False):

        P = np.zeros((0, 2))
        P = np.append(P, self.line_segment(self.Au, self.Ad, N), axis=0)
        P = np.append(P, self.line_segment(self.Bu, self.Bd, N), axis=0)
        P = np.append(P, self.line_segment(self.Au, self.Bu, N), axis=0)
        P = np.append(P, self.line_segment(self.Ad, self.Cl, N), axis=0)
        P = np.append(P, self.line_segment(self.Cr, self.Bd, N), axis=0)
        P = np.append(P, self.line_segment(self.Cl, self.Dl, N), axis=0)
        P = np.append(P, self.line_segment(self.Cr, self.Dr, N), axis=0)
        P = np.append(P, self.line_segment(self.Dl, self.Dr, N), axis=0)

        if(DISPLAY):
            plt.scatter(P.T[0], P.T[1], s=1)
            plt.show()

        return P


    #特徴点群を返す
    def featurePcd(self):

        P = np.zeros((0, 2))
        P = np.append(P, [self.Au], axis=0)
        P = np.append(P, [self.Ad], axis=0)
        P = np.append(P, [self.Bu], axis=0)
        P = np.append(P, [self.Bd], axis=0)
        P = np.append(P, [self.Cl], axis=0)
        P = np.append(P, [self.Cr], axis=0)
        P = np.append(P, [self.Dl], axis=0)
        P = np.append(P, [self.Dr], axis=0)

        return P

    #点群の特徴点を推定する
    def FeaturePointEstimation(self, DISPLAY=False):
        
        if(DISPLAY):
            print("plot 2D...")
            plt.scatter(self.p.T[0], self.p.T[1])
            plt.show()

        self.minimaze()
        self.grouping()
        self.contourEstimation()

        if(DISPLAY):
            self.resultPlot()
            self.resultFeature()


    #T型目的関数
    def Tobj_func(self, p, xy): 
    
        #パラメータ p = [Cx, Cy, α, β] xy = [x, y]
        SS = 0
        
        #直線の係数(ax + by + c = 0)
        tana = np.tan(p[2])
        tanb = np.tan(p[3])
        I = np.array([-1, tana, p[0] - tana*p[1]])
        II = np.array([tanb, -1, p[1] - tanb*p[0]])
        
        
        for i in range(len(xy)):
            
            ss = np.zeros(2)
            #点と直線の距離を求める
            ss[0] = (I[0]*xy[i][0] + I[1]*xy[i][1] + I[2])**2 / (I[0]**2 + I[1]**2)
            ss[1] = (II[0]*xy[i][0] + II[1]*xy[i][1] + II[2])**2 / (II[0]**2 + II[1]**2)
            
            SS += np.min(ss)
            
        return SS


    #フィッティング
    def minimaze(self):

        #T型初期値の設定
        #p = [Cx, Cy, Θ]
        p = np.array([np.mean(self.p.T[0]), np.max(self.p.T[1]), 0, 0])
        P = minimize(self.Tobj_func, p, args=self.p, method="powell")

        # 結果を格納
        tana = np.tan(P.x[2])
        tanb = np.tan(P.x[3])
        self.I = np.array([[-1, tana, P.x[0] - tana*P.x[1]], [tanb, -1, P.x[1] - tanb*P.x[0]]])
        

    #点群をグルーピング
    def grouping(self):

        def distance_cal(a, b, xy):
            d = b - a
            # 線分のベクトル方程式:　a + td のtを求める
            t = np.dot((xy - a), d) / (np.linalg.norm(d)**2)
            
            return np.linalg.norm(xy - (a + d*t))**2

        for i in range(len(self.p)):
            
            #点と直線の距離を求める
            Id = (self.I[0][0]*self.p[i][0] + self.I[0][1]*self.p[i][1] + self.I[0][2])**2 / (self.I[0][0]**2 + self.I[0][1]**2)
            IId = (self.I[1][0]*self.p[i][0] + self.I[1][1]*self.p[i][1] + self.I[1][2])**2 / (self.I[1][0]**2 + self.I[1][1]**2)

            if Id < IId:
                self.Ip = np.append(self.Ip, [self.p[i]], axis=0)
            else:
                self.IIp = np.append(self.IIp, [self.p[i]], axis=0)

        #仮のT端子A, B, C, D点を置く
        self.A = np.array([self.cfg[0][1], -(self.I[1][0]*self.cfg[0][1] + self.I[1][2]) / self.I[1][1]])
        self.B = np.array([self.cfg[0][0], -(self.I[1][0]*self.cfg[0][0] + self.I[1][2]) / self.I[1][1]])
        self.C = np.array([-(self.I[0][1]*self.cfg[1][0] + self.I[0][2]) / self.I[0][0], self.cfg[1][0]])
        self.D = np.array([-(self.I[0][1]*self.cfg[1][1] + self.I[0][2]) / self.I[0][0], self.cfg[1][1]])


        #各点と線分との距離を格納する配列群
        Al, Bl, Cl, Dl = np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0)

        #I(web)要素点群の分類とその距離を求める
        DC = self.C - self.D
        for i in range(len(self.Ip)):

            d = distance_cal(self.D, self.C, self.Ip[i])

            if np.cross(DC, self.Ip[i] - self.D) >= 0:
                self.Cp = np.append(self.Cp, [self.Ip[i]], axis=0)
                Cl = np.append(Cl, [d])
            else:
                self.Dp = np.append(self.Dp, [self.Ip[i]], axis=0)
                Dl = np.append(Dl, [d])

        #II(フランジ)要素点群の分類とその距離を求める
        AB = self.B- self.A
        for i in range(len(self.IIp)):

            d = distance_cal(self.A, self.B, self.IIp[i])

            if np.cross(AB, self.IIp[i] - self.A) >= 0:
                self.Ap = np.append(self.Ap, [self.IIp[i]], axis=0)
                Al = np.append(Al, [d])
            else:
                self.Bp = np.append(self.Bp, [self.IIp[i]], axis=0)
                Bl = np.append(Bl, [d])
        
        #=============================四分位範囲処理=============================
        al = [np.percentile(Al, 25), np.percentile(Al, 75)] #25, 75パーセンタイルを格納
        bl = [np.percentile(Bl, 25), np.percentile(Bl, 75)]
        cl = [np.percentile(Cl, 25), np.percentile(Cl, 75)]
        dl = [np.percentile(Dl, 25), np.percentile(Dl, 75)]

        #四分位範囲を満たす点群を格納する配列群
        ap, bp, cp, dp = np.zeros((0, 2)), np.zeros((0, 2)), np.zeros((0, 2)), np.zeros((0, 2))
            
        for i in range(len(self.Ap)):
            if al[0] <= Al[i] and Al[i] <= al[1]:
                ap = np.append(ap, [self.Ap[i]], axis=0)

        for i in range(len(self.Bp)):
            if bl[0] <= Bl[i] and Bl[i] <= bl[1]:
                bp = np.append(bp, [self.Bp[i]], axis=0)

        for i in range(len(self.Cp)):
            if cl[0] <= Cl[i] and Cl[i] <= cl[1]:
                cp = np.append(cp, [self.Cp[i]], axis=0)
        
        for i in range(len(self.Dp)):
            if dl[0] <= Dl[i] and Dl[i] <= dl[1]:
                dp = np.append(dp, [self.Dp[i]], axis=0)

        #更新
        self.Ap = ap
        self.Bp = bp
        self.Cp = cp
        self.Dp = dp


    def resultPlot(self):

        plt.scatter(self.Ip.T[0], self.Ip.T[1])
        plt.scatter(self.IIp.T[0], self.IIp.T[1])
        plt.show()

        plt.scatter(self.Ap.T[0], self.Ap.T[1])
        plt.scatter(self.Bp.T[0], self.Bp.T[1])
        plt.scatter(self.Cp.T[0], self.Cp.T[1])
        plt.scatter(self.Dp.T[0], self.Dp.T[1])
        plt.show()


    def resultFeature(self):

        plt.scatter(self.Au[0], self.Au[1])
        plt.scatter(self.Ad[0], self.Ad[1])
        plt.scatter(self.Bu[0], self.Bu[1])
        plt.scatter(self.Bd[0], self.Bd[1])
        plt.scatter(self.Cl[0], self.Cl[1])
        plt.scatter(self.Cr[0], self.Cr[1])
        plt.scatter(self.Cl[0], self.Cl[1])
        plt.scatter(self.Dr[0], self.Dr[1])
        plt.scatter(self.Dl[0], self.Dl[1])
        plt.show()
    

    # フランジ、ウェブの輪郭線分を推定し、特徴点を求める
    def contourEstimation(self):
        
        # 輪郭線分係数行列(Abu, ABu, CDl, CDr, AA, BB, DDの順)
        # y = ax + b or x = ay + b --> [a, b]
        I = np.zeros((7, 2)) 

        bsm = BiPlelLineSegModel(self.Ap, self.Bp)
        s = bsm.yminimize()
        I[0] = np.array([s[0], s[1]])
        I[1] = np.array([s[0], s[2]])

        bsm = BiPlelLineSegModel(self.Cp, self.Dp)
        s = bsm.xminimize()
        I[2] = np.array([s[0], s[1]])
        I[3] = np.array([s[0], s[2]])


        #ACの距離をdとすれば、d*5/6以上の距離にある点群をAA線分の要素点群とする（他も同様のアルゴリズムを採用する)
        AC = np.linalg.norm(self.C - self.A)
        BC = np.linalg.norm(self.C - self.B)
        DC = np.linalg.norm(self.C - self.D)

        dc = self.D - self.C

        #各端子点群
        pa, pb, pd = np.zeros((0, 2)), np.zeros((0, 2)), np.zeros((0, 2))

        for i in range(len(self.IIp)):

            d = np.linalg.norm(self.IIp[i] - self.C)
            param = np.cross(dc, self.IIp[i] - self.D) #外積の正負で左右を判別

            if (d > AC*(5/6)) and (param > 0):
                pa = np.append(pa, [self.IIp[i]], axis=0)

            elif (d > BC*(5/6)) and (param < 0):
                pb = np.append(pb, [self.IIp[i]], axis=0)


        for i in range(len(self.Ip)):

            if np.linalg.norm(self.Ip[i] - self.C) > DC*(9/10):
                pd = np.append(pd, [self.Ip[i]], axis=0)
        

        bsm = BiPlelLineSegModel(pa, pb)
        s = bsm.xminimize() #最適化
        I[4] = np.array([s[0], s[1]])
        I[5] = np.array([s[0], s[2]])

        tdm = TwoDimlinearModel(pd)
        s = tdm.linearCoef_1()
        I[6] = np.array([s[0], s[1]])


        #特徴点を求める
        Auy = (I[0][0]*I[4][1] + I[0][1]) / (1 - I[0][0]*I[4][0])
        Ady = (I[1][0]*I[4][1] + I[1][1]) / (1 - I[1][0]*I[4][0])
        Buy = (I[0][0]*I[5][1] + I[0][1]) / (1 - I[0][0]*I[5][0])
        Bdy = (I[1][0]*I[5][1] + I[1][1]) / (1 - I[1][0]*I[5][0])
        Cly = (I[1][0]*I[2][1] + I[1][1]) / (1 - I[1][0]*I[2][0])
        Cry = (I[1][0]*I[3][1] + I[1][1]) / (1 - I[1][0]*I[3][0])
        Dly = (I[6][0]*I[2][1] + I[6][1]) / (1 - I[6][0]*I[2][0])
        Dry = (I[6][0]*I[3][1] + I[6][1]) / (1 - I[6][0]*I[3][0])

        self.Au = np.array([I[4][0]*Auy + I[4][1], Auy])
        self.Ad = np.array([I[4][0]*Ady + I[4][1], Ady])
        self.Bu = np.array([I[5][0]*Buy + I[5][1], Buy])
        self.Bd = np.array([I[5][0]*Bdy + I[5][1], Bdy])
        self.Cl = np.array([I[2][0]*Cly + I[2][1], Cly])
        self.Cr = np.array([I[3][0]*Cry + I[3][1], Cry])
        self.Dl = np.array([I[2][0]*Dly + I[2][1], Dly])
        self.Dr = np.array([I[3][0]*Dry + I[3][1], Dry])


    #始点から終点間をN分割したときの各々の点座標を求める
    def line_segment(self, p1, p2, N): #p1, p2は始点と終点

        P = np.zeros((N, p1.shape[0]))
        
        for i in range(len(P)):
            P[i] = p1 + ((i + 1)/N)*(p2 - p1)
        
        return P


# 二次元線形モデル
class TwoDimlinearModel:

    xy = np.zeros((0, 2))
    Cov = np.zeros((2, 2)) #分散共分散行列
    Avg = np.zeros(2) #平均

    def __init__(self, p):
        self.xy = p
        self.Avg = np.array([np.mean(self.xy.T[0]), np.mean(self.xy.T[1])])
        self.Cov = np.cov(self.xy.T, bias=True)

    #ax - y + b = 0　の回帰係数を返す
    def linearCoef_1(self):

        a = self.Cov[0][1] / self.Cov[0][0]
        b = self.Avg[1] - a*self.Avg[0]

        return a, b

    #ay - x + b = 0 の回帰係数を返す
    def linearCoef_2(self):

        a = self.Cov[0][1] / self.Cov[1][1]
        b = self.Avg[0] - a*self.Avg[1]

        return a, b


# 二平行線分推定モデル
class BiPlelLineSegModel:

    Pa = np.zeros((0, 2))
    Pb = np.zeros((0, 2))

    def __init__(self, Pa, Pb):
        self.Pa = Pa
        self.Pb = Pb


    #最小二乗法による推定(y = ax + b)
    def yminimize(self):

        d, c = np.zeros(2), np.zeros(2)

        tdm = TwoDimlinearModel(self.Pa)
        d[0], c[0] = tdm.linearCoef_1()
        tdm = TwoDimlinearModel(self.Pb)
        d[1], c[1] = tdm.linearCoef_1()

        # パラメータ初期値(d, c1, c2)
        pm = np.array([np.mean(d), c[0], c[1]])
        P = minimize(self.yobj_func, pm, args=[self.Pa, self.Pb], method="powell")

        return P.x[0], P.x[1], P.x[2]

    
    # 最小二乗法による推定(x = ay + b)
    def xminimize(self):

        d, c = np.zeros(2), np.zeros(2)

        tdm = TwoDimlinearModel(self.Pa)
        d[0], c[0] = tdm.linearCoef_2()
        tdm = TwoDimlinearModel(self.Pb)
        d[1], c[1] = tdm.linearCoef_2()

        # パラメータ初期値(d, c1, c2)
        pm = np.array([np.mean(d), c[0], c[1]])
        P = minimize(self.xobj_func, pm, args=[self.Pa, self.Pb], method="powell")

        return P.x[0], P.x[1], P.x[2]


    #目的関数(y = ax + b)
    def yobj_func(self, ps, xyxy):

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


    # 目的関数(x = ya + b)
    def xobj_func(self, ps, xyxy):

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






