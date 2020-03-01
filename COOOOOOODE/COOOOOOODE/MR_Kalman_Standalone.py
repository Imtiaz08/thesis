



from ReadRinexNav import ReadRinexNav 
from ReadRinexRealNav import ReadRinexRealNav
from ReadRinexSprientObs import ReadRinexSprientObs
from ReadRinexRealObs import ReadRinexRealObs
from ComputeOrbit import ComputeOrbit
import math
import numpy as np
from numpy import mat
from scipy.linalg import block_diag
import os
import pandas
from datetime import datetime as dt
import matplotlib.pyplot as plt
from Coordinateconversion import*
import pdb
from collections import Counter

def Kalman_Standalone(obs, ephemeris, car):

    XYZ_Station = obs[-1]
    
    
    # Set Q state transition variance
#    Sf = 0.01; Sg = 0.01; sigma=0.01;                      
    Sf = 36;             # 36 360 3600 36000 360000
    Sg = 0.01; sigma=1;                      
    Qb = np.mat([[Sf  + Sg / 3, Sg / 2], [Sg / 2, Sg]])                            
    Qxyz = sigma ** 2 * np.mat([[1 / 3, 1 / 2, 1, 1], [1 / 2, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]])
    Q = block_diag(Qxyz,Qxyz,Qxyz,Qb)                                                       
    
    epochs = [j[1] for i, j in enumerate(obs[0:-1])]
    epoch = list(set(epochs))
    epoch.sort(key=epochs.index)
    r = []
    clock = []
    GDOP_list = []
    WGDOP_list = []
    EDOP_list = []
    NDOP_list = []
    NLOS = {}
    Zglobal = {}
    Zlocal = {}
    test_err = []
    
    for t, tow in enumerate(epoch):
        if t == 0:
            # Set initial values of X and P     
            pos = np.zeros(shape = (14,1))
            
            pos[[0, 4, 8]] = [[XYZ_Station[0]], [XYZ_Station[1]], [XYZ_Station[2]]]               #Initial position
            pos[[1, 5, 9]] = [[0], [0], [0]]
            pos[[2, 6, 10]] = [[0], [0], [0]]
            pos[[3, 7, 11]] = [[0], [0], [0]]                 #Initial velocity
            pos[12] = -1000000            #1200000                              #Clock Bias
            pos[13] = 150                                          #Clock Drift
            P = np.eye(14)*1;
            r.append([])
            clock.append([])
            r[t].append(tow)
            r[t].append(float(pos[0]))
            r[t].append(float(pos[4]))
            r[t].append(float(pos[8]))
            clock[t].append(float(pos[12]))
            clock[t].append(float(pos[13]))
            continue
        
        indx1 = [i for i, j in enumerate(obs) if j[1] == tow]              # indx1: obs 
#         if the car is moving
        indx2 = [i for i, j in enumerate(car) if j[1] == tow]              # indx2: car 
        if indx1 == [] or indx2 == []:
            continue
        pos[1:10:4] = [[float(car[indx2[0]][5])], [float(car[indx2[0]][6])], [float(car[indx2[0]][7])]]             # apply velocity from car umt
        pos[2:11:4] = [[float(car[indx2[0]][8])], [float(car[indx2[0]][9])], [float(car[indx2[0]][10])]]            # apply acceleration from car umt
        pos[3:12:4] = [[0], [0], [0]]                                                                        
#        # if the car is static
#        pos[1:10:4] = [[0], [0], [0]]
#        pos[2:11:4] = [[0], [0], [0]]
#        pos[3:12:4] = [[0], [0], [0]]
  
#==============================================================================
#       print(t, tow)
#        if tow == 325367:
#            pdb.set_trace()
#            print(t)
#            print([int(np.mod(tow, 3600*24)/3600), int(np.mod(np.mod(tow, 3600*24), 3600)/60), np.mod(np.mod(np.mod(tow, 3600*24), 3600),60)])
              
# 
#==============================================================================
        
        Satclocrr = []; X = []; Y = []; Z = []; PR = []; azim = []; elev = []; PRN = [];CNR = []
        
        for i in indx1:
            satorbit = ComputeOrbit(obs[i],ephemeris[0:-2],XYZ_Station)
            if obs[i][6] >= 20 and obs[i][3] > 0:  #Select data with high S/N
                azimuth, elevation = efix2topo(satorbit[4], satorbit[5], satorbit[6], XYZ_Station)
                if elevation >= 5:
                    Satclocrr.append(satorbit[3])
                    PRN.append(satorbit[2])
                    X.append(satorbit[4])
                    Y.append(satorbit[5])
                    Z.append(satorbit[6])
                    CNR.append(obs[i][6])
                    

                    # atmosphere correction
                    
                    T = 2.3 / math.cos((90 - elevation) * math.pi / 180)
                    temp_r = [car[indx2[0]][2], car[indx2[0]][3], car[indx2[0]][4]]
                    I = Ionosphereklobuchar(temp_r,elevation,azimuth,tow,ephemeris[-2],ephemeris[-1])
                    
                    PR.append(obs[i][3] + satorbit[3]- T - I)
                    azim.append(azimuth)
                    elev.append(elevation)
                    #sigma = 293 * 0.1 * math.sqrt(0.5) * pow(10, -(obs[i][6]-2*abs(obs[i][6]-53.99))/20)
                    #Rhoerror.append(sigma ** 2)
#                    Rhoerror.append(0.001)
  
#                    
        SNR, PR, X, Y, Z, elev, azim, grp1, PRN = SNRsort(CNR, PR, X, Y, Z, elev, azim, PRN)            
                    
        # Set R  variance of measurement error(pseudorange error)
                    
        #R = np.eye(len(X)) * Rhoerror                     
# =============================================================================
        #positioning using Kalman Filter
                    
        obs_time = [int(np.mod(tow, 3600*24)/3600), int(np.mod(np.mod(tow, 3600*24), 3600)/60), np.mod(np.mod(np.mod(tow, 3600*24), 3600),60)]
        obs_time_str = ' '.join(str(e) for e in obs_time)
        
        
        
#        print(t, tow, len(X))
        
        if len(X) > 3:
            [pos,P,DOP,WDOP, NLOS_t, Zglobal_t, Zlocal_t, err] = Extended_KF(Q, P, PR, pos, X, Y, Z, grp1, PRN);
            GDOP = math.sqrt(DOP[0,0]+DOP[1,1]+DOP[2,2]+DOP[3,3])
            WGDOP = math.sqrt(WDOP[0,0]+WDOP[1,1]+WDOP[2,2]+WDOP[3,3])
            DOPenu = DOPxyz2enu([car[indx2[0]][2], car[indx2[0]][3], car[indx2[0]][4]] , DOP)
            EDOP_list.append(math.sqrt(DOPenu[0,0]))
            NDOP_list.append(math.sqrt(DOPenu[1,1]))
        else:
            pos, _ = CarVelocity(pos ,1)
            DOP = np.nan
            WDOP = np.nan
            NLOS_t = np.nan
            Zglobal_t = np.nan
            Zlocal_t = np.nan
            err = np.nan
            WGDOP = np.nan
            EDOP_list.append(np.nan)
            NDOP_list.append(np.nan)
            
        #NLOS[(t, tow, obs_time_str)] = NLOS_t
        NLOS[t] = tow, obs_time_str, NLOS_t
        Zglobal[t] = tow, obs_time_str, Zglobal_t
        Zlocal[t] = tow, obs_time_str, Zlocal_t
        test_err.append(err)
# =============================================================================

# =============================================================================
        r.append([])
        clock.append([])
#        print(indx2)
        r[t].append(tow)
        r[t].append(float(pos[0]))
        r[t].append(float(pos[4]))
        r[t].append(float(pos[8]))
        clock[t].append(float(pos[12]))
        clock[t].append(float(pos[13]))
        GDOP_list.append(GDOP)
        WGDOP_list.append(WGDOP)
        XYZ_Station = [float(pos[0]), float(pos[4]), float(pos[8])]
        
    return r, GDOP_list, WGDOP_list, EDOP_list, NDOP_list, clock, NLOS, Zglobal, Zlocal, test_err
        
        
def Extended_KF(Q, P, PR, pos, X, Y, Z, grp1, PRN):
    Ii = np.eye(len(pos))
    Ii = mat(Ii)

    
    Xp, fy = CarVelocity(pos ,1)
    
       
    
    fy = mat(fy)
    P = mat(P)
    
    Pp = fy * P * fy.T + Q
    
    [Gv, Ddb] = Gv_solution([Xp[0], Xp[4], Xp[8], Xp[12]], X, Y, Z)
    
    R, nlos_t, Zg_t, Zl_t, deleteindx = NLOSdect(Gv, Ddb, PR, grp1, PRN)    # R of satellites
    
#==============================================================================
#    if deleteindx and len(X)-len(deleteindx) >= 4:
#        X = np.delete(X, deleteindx) 
#        Y = np.delete(Y, deleteindx)
#        Z = np.delete(Z, deleteindx)
#        PR = np.delete(PR, deleteindx)
#        R = np.delete(R, deleteindx, axis = 0)
#        R = np.delete(R, deleteindx, axis = 1)
#==============================================================================
    
    
    gXp, H = ObservationM(Xp, X, Y, Z)
    
    PR = mat(PR)
    gXp = mat(gXp)
    
    Yk = PR.T - gXp  ####################################################140m 
    
    Pp = mat(Pp)
    H = mat(H)
    R = mat(R)
    
    K = Pp * H.T * (H * Pp * H.T + R).I
    
    Xp = mat(Xp)
    
    Xo = Xp + K * Yk
    
    
    Po = (Ii - K * H) * Pp
    
    # Calculate DOP
    G = H.T[~(H.T == 0).all(1).A1].T
    
    DOP = (G.T * G).I
    
    
    W = np.eye(len(X)) * 1/np.diag(R)
    
    WDOP = (G.T * W * G).I

    
    return Xo, Po, DOP, WDOP, nlos_t, Zg_t, Zl_t, np.mean(Yk)

        
def CarVelocity(X, T):
    Val = np.zeros(shape = (14,1))
    Val[0:9:4] = X[0:9:4] +  T * X[1:10:4] + 0.5 * T ** 2 * X[2:11:4] + T ** 3 * X[3:12:4]/6    # POS
    Val[1:10:4] = X[1:10:4] + T * X[2:11:4] + 0.5 * T ** 2 * X[3:12:4]                          # VELOCITY
    Val[2:11:4] = X[2:11:4] + T * X[3:12:4]                                                     # ACCELERATION
    Val[3:12:4] = X[3:12:4]                                                                     # jerk 
    Val[12] = X[12] + T * X[13]                                                                 #clock bias
    Val[13] = X[13]                                                                             # drift
    Jacob1 = np.mat([[1, T, 0.5 * T ** 2, T ** 3 / 6], [0, 1, T, 0.5 * T ** 2], [0, 0, 1, T], [0, 0, 0, 1]])
    Jacob2 = np.mat([[1, T], [0, 1]])
    Jacob = block_diag(Jacob1,Jacob1,Jacob1,Jacob2)
    
    return Val, Jacob


def ObservationM(pos, X, Y, Z):
    Val = np.zeros(shape = (len(X), 1))
    dist = np.zeros(shape = (len(X), 1))
    Jacob = np.zeros(shape = (len(X), len(pos)))
    
    for i in range(0, len(X)):
        Val[i] = math.sqrt((X[i] - pos[0]) ** 2 + (Y[i] - pos[4]) ** 2 + (Z[i] - pos[8]) ** 2)
        dist[i] = Val[i] + pos[12]
#        dist[i] = Val[i]
        Jacob[i][0] = (pos[0] - X[i])/Val[i]
        Jacob[i][4] = (pos[4] - Y[i])/Val[i]
        Jacob[i][8] = (pos[8] - Z[i])/Val[i]
        Jacob[i][12] = 1
    return dist, Jacob




def SNRsort(SNR, corr_meas, X, Y, Z, ele, azim, SVID):
    thre = 42 

    SNR_phi1 = []
    corr_meas_phi1 = []
    X_phi1 = []
    Y_phi1 = []
    Z_phi1 = []
    ele_phi1 = []
    azim_phi1 = []
    SVID_phi1 = []

    SNR_phi2 = []
    corr_meas_phi2 = []
    X_phi2 = []
    Y_phi2 = []
    Z_phi2 = []
    ele_phi2 = []
    azim_phi2 = []                
    SVID_phi2 = []    
    
    cnt = 0

    for j, k in enumerate(SNR):
        if k < thre:
            SNR_phi1.append(SNR[j])
            corr_meas_phi1.append(corr_meas[j])
            X_phi1.append(X[j])
            Y_phi1.append(Y[j])
            Z_phi1.append(Z[j])
            ele_phi1.append(ele[j])
            azim_phi1.append(azim[j])
            SVID_phi1.append(SVID[j])            
            
            cnt += 1
        else:
            
            SNR_phi2.append(SNR[j])
            corr_meas_phi2.append(corr_meas[j])
            X_phi2.append(X[j])
            Y_phi2.append(Y[j])
            Z_phi2.append(Z[j])
            ele_phi2.append(ele[j])
            azim_phi2.append(azim[j])                       
            SVID_phi2.append(SVID[j])            
                
    SNR_sorted = np.hstack((SNR_phi1,SNR_phi2))
    corr_meas_sorted = np.hstack((corr_meas_phi1,corr_meas_phi2))
    X_sorted = np.hstack((X_phi1,X_phi2))
    Y_sorted = np.hstack((Y_phi1,Y_phi2))
    Z_sorted = np.hstack((Z_phi1,Z_phi2))
    ele_sorted = np.hstack((ele_phi1,ele_phi2))
    azim_sorted = np.hstack((azim_phi1,azim_phi2))
    SVID_sorted = np.hstack((SVID_phi1,SVID_phi2))
    
    return SNR_sorted, corr_meas_sorted, X_sorted, Y_sorted, Z_sorted, ele_sorted, azim_sorted, cnt, SVID_sorted
    
    
   
def NLOSdect(Gv, D, meas, L1, PRN):
#==============================================================================
#     [nrows, ncols] = np.shape(H)
#     Gv = np.zeros(shape = (nrows, 1))
#     for i in range(0, num_car):       
#         Gv = np.hstack((Gv, H[:,14*i:14*i+12:4]))
#     
#     Gv = np.delete(Gv, 0, axis = 1)           # 13*6
#==============================================================================
#    sigma = 293 * 0.1 * math.sqrt(0.5) * pow(10, -(snr-2*abs(snr-53.99))/20)
#    Rhoerrorp = sigma ** 2
#    
    num_sate = len(meas)
    Rhoerrorp = 9                                           
    R = np.eye(num_sate) * Rhoerrorp                          # measurement noise of PR
       
        
    delta = np.eye(4)                              # predicted position error #########################333
    r_nlos = 10
    
    kesi_V = np.matmul(np.matmul(Gv, delta), Gv.T) + R            # squared matrix
    kesi_V = mat(kesi_V)
    V11 = kesi_V[0:L1, 0:L1]
    V21 = kesi_V[L1:, 0:L1]
    V12 = kesi_V[0:L1, L1:]
    V22 = kesi_V[L1:, L1:]
    
    
    D1 = D[:L1] #+ np.ones((L1,1)) * r_nlos        ##############################################
    D2 = D[L1:]
    D1 = mat(D1)
    D2 = mat(D2)
    
    meas_L1 = meas[:L1]
    meas_L2 = meas[L1:]
    meas_L2 = mat(meas_L2)
    meas_L1 = mat(meas_L1)
    
    #miu = D1 + V12 * np.linalg.inv(V22) * (meas_L2 - D2)
    miu = D1 + np.matmul(np.matmul(V12, np.linalg.inv(V22)), (meas_L2.T - D2))  # L1 * 1         mostly D1 
    
    cov_v = V11 - np.matmul(np.matmul(V12, np.linalg.inv(V22)), V21)   
    
    Z_zeta = (meas_L1.T - miu).T * np.linalg.inv(cov_v) * (meas_L1.T - miu)   # 1*1
       
    
    Z_itrt = {}
    NLOS_sig = {}
    sate_deleteindx = []    
    
    if L1 > 0 and Z_zeta > 5:      # global threshold 
        for i in range(L1):
            D1 = D[i] #+ r_nlos * np.ones((1,1))    ################## ????????????????
            D2 = D[L1:]
            
            meas_L1 = meas[i]
            meas_L2 = meas[L1:]
            meas_L2 = mat(meas_L2)
            meas_L1 = mat(meas_L1)
            
            V11 = kesi_V[i,i]
            V12 = kesi_V[i, L1:]
            V21 = kesi_V[L1:, i]
                
            miu = D1 + np.matmul(np.matmul(V12, np.linalg.inv(V22)), (meas_L2.T - D2))
    
            cov_v = V11 - np.matmul(np.matmul(V12, np.linalg.inv(V22)), V21)   
            
            Z = (meas_L1.T - miu).T * np.linalg.inv(cov_v) * (meas_L1.T - miu)
           
            if Z >= 2:
                R[i, i] = Z * R[i, i] 
                NLOS_sig[PRN[i]] = meas_L1
                Z_itrt[PRN[i]] = Z
                sate_deleteindx.append(i)
                        
    RR = R[:num_sate, :num_sate]   # only extract sate measurements
    
    return RR, NLOS_sig, Z_zeta, Z_itrt, sate_deleteindx 

    

  
    
def Gv_solution(pre_pos,X,Y,Z):
    
    dist = np.zeros(shape = (len(X), 1))
    G_sate = np.zeros(shape = (len(X), 4))          #############################################3
    dist_prime = np.zeros(shape = (len(X), 1))      # true range plus clk bias
    
    
    for i in range(0, len(X)):
        dist[i] = np.sqrt((X[i] - pre_pos[0]) ** 2 + (Y[i] - pre_pos[1]) ** 2 + (Z[i] - pre_pos[2]) ** 2)     # sqrt(x^2+y^2+z^2) denominator
        dist_prime[i] = dist[i] + pre_pos[3]            # clk bias 
        G_sate[i][0] = (pre_pos[0] - X[i])/dist[i]
        G_sate[i][1] = (pre_pos[1] - Y[i])/dist[i]
        G_sate[i][2] = (pre_pos[2] - Z[i])/dist[i]
        G_sate[i][3] = 1                             ####################################
        
    
    return G_sate, dist_prime



def Date2GPSTime(utc):
    gps_week_start = dt(1980,1,6,0,0,0)
    date_time = dt(int(utc[0]),int(utc[1]),int(utc[2]),int(utc[3]),int(utc[4]),int(float(utc[5])))
    
    gpsstartdt = 366 + gps_week_start.toordinal() + (gps_week_start - dt.fromordinal(gps_week_start.toordinal())).total_seconds()/(24*60*60)
    utcdt = 366 + date_time.toordinal() + (date_time - dt.fromordinal(date_time.toordinal())).total_seconds()/(24*60*60)
    
    tmp = (utcdt - gpsstartdt)/7
    GPS_Weeks = math.floor(tmp)
    GPS_SOW = round((tmp-GPS_Weeks)*7*24*3600)
    
    return GPS_Weeks, GPS_SOW



def ExtractCar(car, flag):
    if flag == 1:
        car = np.array(car)
        car = car.tolist()
        year = 2019;
        month = 7;
        day = 3;
        extract_car_all = []
        for i in range(0,len(car)):
            if int(car[i][0]) - car[i][0] == 0:
                extract_car = []
                hour = int(car[i][0]/3600)
                minute = int((car[i][0]-3600*hour)/60)
                second = car[i][0]-3600*hour-60*minute
                utc = [year,month,day,hour,minute,second]
                GPS_Weeks, GPS_SOW = Date2GPSTime(utc)
                GPS_SOW = GPS_SOW + 6*3600
                extract_car.append(GPS_Weeks)
                extract_car.append(GPS_SOW)
                extract_car.append(car[i][3])
                extract_car.append(car[i][4])
                extract_car.append(car[i][5])
                extract_car.append(car[i][6])
                extract_car.append(car[i][7])
                extract_car.append(car[i][8])
                extract_car.append(car[i][9])
                extract_car.append(car[i][10])
                extract_car.append(car[i][11])
                extract_car_all.append(extract_car)
                
    elif flag == 2:
        car = np.array(car)
        car = car.tolist()
        extract_car_all = []
        for i in range(0,len(car)):
            if int(car[i][0]) - car[i][0] == 0:
                extract_car = []
                extract_car.append(2000) # whatever now
                extract_car.append(car[i][0])
                extract_car.append(car[i][3])
                extract_car.append(car[i][4])
                extract_car.append(car[i][5])
                extract_car.append(car[i][6])
                extract_car.append(car[i][7])
                extract_car.append(car[i][8])
                extract_car.append(car[i][9])
                extract_car.append(car[i][10])
                extract_car.append(car[i][11])
                extract_car_all.append(extract_car)
                    
    return extract_car_all






if __name__ == '__main__':
#    crf = os.getcwd()
#    ephemeris, _, _, _, _, = ReadRinexNav(crf+'/Data/Nav/03072019.nav')
#    gps, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0000.obs')
#    car = pandas.read_csv(crf + '/Data/umt/v2.umt', header = None); car = ExtractCar(car)
    
#    gps, _, _, _ = ReadRinexRealObs("./meas/Obs/190523_1.obs")
#    ephemeris, _, _, _, _, = ReadRinexRealNav("./meas/Nav/190523_1.nav")
#    car = pandas.read_csv("./meas/Car/190523_1.umt", header = None); car = ExtractCar(car,2)
    
    
    gps, _, _, _ = ReadRinexRealObs("./meas/Obs/newVeh0481.obs")
    ephemeris, _, _, _, _, = ReadRinexRealNav("./meas/Nav/03072019.nav")
    car = pandas.read_csv("./meas/Car/newVeh0481.umt", header = None); car = ExtractCar(car,1)    
    
    
    
    r, GDOP, WGDOP, EDOP, NDOP, clock, NLOSDection, statsglobal, statslocal, test_err = Kalman_Standalone(gps,ephemeris,car)
    
    print(np.nanmean(test_err))
    gg = list(list(zip(* list(statsglobal.values())))[2])
    gg = np.vstack(gg[1:])                     #################### the first one is excluded
    gg_NLOS_cnt = sum(i > 5 for i in gg)   
    print(gg_NLOS_cnt)
    plt.figure()
    plt.plot(gg[250:,0], label = 'Global test')
    plt.legend(loc = 'lower right')
    
    
    
    est_clk = np.vstack(list(zip(* clock))[0])    
    clk_dft = np.vstack(list(zip(* clock))[1])
    plt.figure()
    plt.plot(est_clk, label = 'Receiver Clock offset')
    plt.plot(clk_dft, label = 'Receiver Clock drift', color = 'r')
    plt.legend(loc='lower right')
    plt.xlabel('Time(sec)')
    plt.ylabel('m')    
    
    
#    est_car_x = mat(list(zip(* r))[1])
#    est_car_y = mat(list(zip(* r))[2])
#    est_car_z = mat(list(zip(* r))[3])    
#
#    ego_car_x = mat(list(zip(* car))[2])
#    ego_car_y = mat(list(zip(* car))[3])
#    ego_car_z = mat(list(zip(* car))[4])    
    
    
    
    time = list(set(list(zip(* car))[1]) & set(list(zip(* r))[0]))    
    
    startpos = list(list(zip(* car))[1]).index(time[0])
    car_origin = [car[startpos][2],car[startpos][3], car[startpos][4]]       #### static case
    lon0, lat0, h0 = Geodetic(car_origin) 
    
    est_car = [[val[1], val[2], val[3]] for i, val in enumerate(r) if val[0] in time ]
    ego_car = [[val[2], val[3], val[4]] for i, val in enumerate(car) if val[1] in time ]
    est_car = mat(est_car)
    ego_car = mat(ego_car)
    
    xEast1, yNorth1, zUp1 = ecef_to_enu(ego_car[:,0], ego_car[:,1], ego_car[:,2], lat0, lon0, h0)        
    xEast2, yNorth2, zUp2 = ecef_to_enu(est_car[:,0], est_car[:,1], est_car[:,2], lat0, lon0, h0)    
    
    
    diff_prime = [xEast2 - xEast1, yNorth2 - yNorth1, zUp2 - zUp1]
    rmse = np.zeros(3)    
    rmse[0] = RMSE(xEast2, xEast1)
    rmse[1] = RMSE(yNorth2, yNorth1)    
    rmse[2] = RMSE(zUp2, zUp1)  
    print(rmse)
   
   
   
#    t = range(np.size(diff_prime[0]))   
          
    fig, axs = plt.subplots(3, sharex = True)
    fig.suptitle('ENU deviation')
    axs[0].plot(diff_prime[0][250:,0], color = 'r', label = 'xEast')
    axs[1].plot(diff_prime[1][250:,0], color = 'r', label = 'yNorth')
    axs[2].plot(diff_prime[2][250:,0], color = 'r', label = 'zUp')
#    axs[0].set_ylim([-5, 10])
#    axs[1].set_ylim([-10, 20])
#    axs[2].set_ylim([-20, 60])
    axs[0].legend()
    axs[1].legend()    
    axs[2].legend()
    plt.xlabel('Time(sec)')
    plt.ylabel('m')
    
    plt.show()
#==============================================================================
#     
#     diff = [est_car_x - ego_car_x, est_car_y - ego_car_y, est_car_z - ego_car_z]
#           
#     plt.figure()
#     plt.plot(diff[0].T, label = 'x')
#     plt.plot(diff[1].T, label = 'y')
#     plt.plot(diff[2].T, label = 'z')
#     plt.legend(loc='lower right')
#     plt.xlabel('Time(sec)')
#     plt.ylabel('m')
#==============================================================================
        

    #plt.close('all')
#    file = open('spp.txt', 'w')
#    for num, line in enumerate(r):
#        file.write(str(r[num][0]))
#        file.write('\t')
#        file.write(str(r[num][1]))
#        file.write('\t')
#        file.write(str(r[num][2]))
#        file.write('\t')
#        file.write(str(r[num][3]))
#        file.write('\n')
#    file.close()
#    
#    file1 = open('GDOP.txt', 'w')
#    for num, line in enumerate(GDOP):
#        file1.write(str(GDOP[num]))
#        file1.write('\n')
#    file1.close()
#    
#    file1 = open('WGDOP.txt', 'w')
#    for num, line in enumerate(WGDOP):
#        file1.write(str(WGDOP[num]))
#        file1.write('\n')
#    file1.close()
#    
#    file1 = open('EDOP.txt', 'w')
#    for num, line in enumerate(EDOP):
#        file1.write(str(EDOP[num]))
#        file1.write('\n')
#    file1.close()
#    
#    file1 = open('NDOP.txt', 'w')
#    for num, line in enumerate(NDOP):
#        file1.write(str(NDOP[num]))
#        file1.write('\n')
#    file1.close()
#    
#    file1 = open('Clock.txt', 'w')
#    for num, line in enumerate(clock):
#        file1.write(str(clock[num][0]))
#        file1.write('\n')
#    file1.close()