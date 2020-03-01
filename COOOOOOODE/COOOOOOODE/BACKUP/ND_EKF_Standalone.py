



from ReadRinexNav import ReadRinexNav 
from ReadRinexSprientObs import ReadRinexSprientObs
from ComputeOrbit import ComputeOrbit
import math
import numpy as np
from numpy import mat
from scipy.linalg import block_diag

def Kalman_Standalone(obs, ephemeris):

    XYZ_Station = obs[-1]
    T = 1 # positioning interval
    
    
    # Set Q state transition variance
    Sf = 36.0; Sg = 1; sigma=5.0;                               # ????????????????????
    Qb = np.mat([[Sf * T + Sg * T ** 3 / 3, Sg * T ** 2 / 2], [Sg * T ** 2 / 2, Sg * T]])                            
    Qxyz = sigma ** 2 * np.mat([[T ** 3 / 3, T ** 2 / 2, T, 1], [T ** 2 / 2, T, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]])
    Q = block_diag(Qxyz,Qxyz,Qxyz,Qb)
    
    epochs = [j[1] for i, j in enumerate(obs[0:-1])]
    epoch = list(set(epochs))
    epoch.sort(key=epochs.index)
    r = []
    NLOS = {}
    for t, tow in enumerate(epoch):
        if t == 0:
            # Set initial values of X and P     
            pos = np.zeros(shape = (14,1))
            
            pos[[0, 4, 8]] = [[XYZ_Station[0]], [XYZ_Station[1]], [XYZ_Station[2]]]               #Initial position
            pos[[1, 5, 9]] = [[0], [0], [0]]
            pos[[2, 6, 10]] = [[0], [0], [0]]
            pos[[3, 7, 11]] = [[0], [0], [0]]                 #Initial velocity
            pos[12] = 0                                       #Clock Bias
            pos[13] = 0                                       #Clock Drift
            P = np.eye(14)*1;                                 
            r.append([])
            r[t].append(float(pos[0]))
            r[t].append(float(pos[4]))
            r[t].append(float(pos[8]))
            continue
        
        
        indx1 = [i for i, j in enumerate(obs[0:-1]) if j[1] == tow]    # satellites
        
        ################## if the car is moving ######################
#        indx2 = [i for i, j in enumerate(obs[0:-1]) if j[1] == tow]
#        if indx1 == [] or indx2 == []:
#            continue
#        pos[1:10:4] = [[float(car[indx2,5])], [float(car[indx2,6])], [float(car[indx2,7])]]                   
#        pos[2:11:4] = [[float(car[indx2,8])], [float(car[indx2,9])], [float(car[indx2,10])]]                  
#        pos[3:12:4] = [[float(car[indx2,11])], [float(car[indx2,12])], [float(car[indx2,13])]]
        
        
        
        ################## if the car is static #########################
        pos[1:10:4] = [[0], [0], [0]]           # velocity
        pos[2:11:4] = [[0], [0], [0]]           # acceleration
        pos[3:12:4] = [[0], [0], [0]]           # jerk
        
        Satclocrr = []; X = []; Y = []; Z = []; PR = []; azim = []; elev = [];snr = []
        for i in indx1:
            satorbit = ComputeOrbit(obs[i],ephemeris[0:-2],XYZ_Station)
            if obs[i][6] >= 0:  #Select data with high S/N
                azimuth, elevation = efix2topo(satorbit[4], satorbit[5], satorbit[6], XYZ_Station)
                if elevation >= 10:
                    Satclocrr.append(satorbit[3])
                    X.append(satorbit[4])
                    Y.append(satorbit[5])
                    Z.append(satorbit[6])
                    PR.append(obs[i][3] + satorbit[3])
                    azim.append(azimuth)
                    elev.append(elevation)
                    snr.append(obs[i][6])
        
        
        
        
        ######################################################################        
        
        R, Z_dect_coupled, Z_dect_decpl, NLOS_detect, SNR_sort, PR, X, Y, Z, ele_sort, azim_sort =\
        SNRdect(snr, PR, X, Y, Z, elev, azim, pos)

        NLOS[t] = NLOS_detect       
        #positioning using Kalman Filter
        [pos,P] = Extended_KF(Q, R, P, PR, pos, X, Y, Z);
        
        r.append([])
        r[t].append(float(pos[0]))
        r[t].append(float(pos[4]))
        r[t].append(float(pos[8]))
        XYZ_Station = [float(pos[0]), float(pos[4]), float(pos[8])]
    return r, NLOS
        
        
def Extended_KF(Q, R, P, PR, pos, X, Y, Z):
    Ii = np.eye(len(pos))
    Ii = mat(Ii)

    
    Xp, _ = CarVelocity(pos ,1)
    
    _, fy = CarVelocity(Xp ,1)
    
    fy = mat(fy)
    P = mat(P)
    
    Pp = fy * P * fy.T + Q
    
    gXp, H = ObservationM(Xp, X, Y, Z)
    
    PR = mat(PR)
    gXp = mat(gXp)
    
    Y = PR.T - gXp
    
    Pp = mat(Pp)
    H = mat(H)
    R = mat(R)
    
    K = Pp * H.T * (H * Pp * H.T + R).I
    
    Xp = mat(Xp)
    
    Xo = Xp + K * Y
    
    
    Po = (Ii - K * H) * Pp
    
    return Xo, Po


def SNRdect(SNR, PR, X, Y, Z, ele, azim, pos):
    thre = 55 
    
    SNR_phi1 = []
    corr_meas_phi1 = []
    X_phi1 = []
    Y_phi1 = []
    Z_phi1 = []
    ele_phi1 = []
    azim_phi1 = []

    SNR_phi2 = []
    corr_meas_phi2 = []
    X_phi2 = []
    Y_phi2 = []
    Z_phi2 = []
    ele_phi2 = []
    azim_phi2 = []                
    
    cnt = 0

    for j, k in enumerate(SNR):
        if k < thre:
            SNR_phi1.append(SNR[j])
            corr_meas_phi1.append(PR[j])
            X_phi1.append(X[j])
            Y_phi1.append(Y[j])
            Z_phi1.append(Z[j])
            ele_phi1.append(ele[j])
            azim_phi1.append(azim[j])
            
            cnt += 1
        else:
            
            SNR_phi2.append(SNR[j])
            corr_meas_phi2.append(PR[j])
            X_phi2.append(X[j])
            Y_phi2.append(Y[j])
            Z_phi2.append(Z[j])
            ele_phi2.append(ele[j])
            azim_phi2.append(azim[j])                       
            
                
    SNR_sorted = np.hstack((SNR_phi1,SNR_phi2))
    corr_meas_sorted = np.hstack((corr_meas_phi1,corr_meas_phi2))
    X_sorted = np.hstack((X_phi1,X_phi2))
    Y_sorted = np.hstack((Y_phi1,Y_phi2))
    Z_sorted = np.hstack((Z_phi1,Z_phi2))
    ele_sorted = np.hstack((ele_phi1,ele_phi2))
    azim_sorted = np.hstack((azim_phi1,azim_phi2))
    

    
    [Gv, Ddb] = Gv_solution(pos, X_sorted,Y_sorted,Z_sorted, corr_meas_sorted)
    
    nb_current_sate = np.size(X)
    
    Rweight, Z, Z1, NLOS_dct = NLOSdect(Gv, Ddb, corr_meas_sorted, cnt, nb_current_sate)    # R of satellites
    

    return Rweight, Z, Z1, NLOS_dct, SNR_sorted, corr_meas_sorted, X_sorted,Y_sorted, Z_sorted,ele_sorted,azim_sorted 
    
    
    
def NLOSdect(Gv, D, meas, L1, num_sate):
#==============================================================================
#     [nrows, ncols] = np.shape(H)
#     Gv = np.zeros(shape = (nrows, 1))
#     for i in range(0, num_car):       
#         Gv = np.hstack((Gv, H[:,14*i:14*i+12:4]))
#     
#     Gv = np.delete(Gv, 0, axis = 1)           # 13*6
#==============================================================================
    Rhoerrorp = 9                                           
    Rp = np.eye(num_sate) * Rhoerrorp                          # measurement noise of PR
       
        
    delta = np.eye(3)
    r_nlos = 10
    
    kesi_V = np.matmul(np.matmul(Gv, delta), Gv.T) + Rp            # squared matrix
    kesi_V = mat(kesi_V)
    V11 = kesi_V[0:L1, 0:L1]
    V21 = kesi_V[L1:, 0:L1]
    V12 = kesi_V[0:L1, L1:]
    V22 = kesi_V[L1:, L1:]
    
    
    D1 = D[:L1] + np.ones((L1,1)) * r_nlos
    D2 = D[L1:]
    
    meas_L1 = meas[:L1]
    meas_L2 = meas[L1:]
    meas_L2 = mat(meas_L2)
    meas_L1 = mat(meas_L1)
    
    #miu = D1 + V12 * np.linalg.inv(V22) * (meas_L2 - D2)
    miu = D1 + np.matmul(np.matmul(V12, np.linalg.inv(V22)), (meas_L2.T - D2))  # L1 * 1
    
    cov_v = V11 - np.matmul(np.matmul(V12, np.linalg.inv(V22)), V21)   
    
    Z_1 = (meas_L1.T - miu).T * np.linalg.inv(cov_v) * (meas_L1.T - miu)   # 1*1
    
    Z_single = []
    NLOS_sig = {}    
    
    thre_Z = 2
    if L1 > 1:
        for i in range(L1):
            D1 = D[i]
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
           
            if Z >= thre_Z:
                Rp[i, i] = Z * Rp[i, i] 
                NLOS_sig[i] = meas_L1
                
            Z_single.append(Z)
                        
    
    return Rp,  Z_1, Z_single, NLOS_sig

    
    
    
    
    
    
def Gv_solution(pos,X,Y,Z,PR_sorted):
    
    dist = np.zeros(shape = (len(X), 1))
    G_sate = np.zeros(shape = (len(X), 3))
    dist_prime = np.zeros(shape = (len(X), 1))      # true range plus clk bias
        
    
    for i in range(0, len(X)):
        dist[i] = np.sqrt((X[i] - pos[0]) ** 2 + (Y[i] - pos[4]) ** 2 + (Z[i] - pos[8]) ** 2)   # sqrt(x^2+y^2+z^2) denominator
        dist_prime[i] = dist[i] + pos[12]
        G_sate[i][0] = (pos[0] - X[i])/dist[i]
        G_sate[i][1] = (pos[4] - Y[i])/dist[i]
        G_sate[i][2] = (pos[8] - Z[i])/dist[i]
    
            
       # dist_prime = dist_prime.tolist()
    
    return G_sate, dist_prime




        
def CarVelocity(X, T):
    Val = np.zeros(shape = (14,1))
    Val[0:9:4] = X[0:9:4] +  T * X[1:10:4] + 0.5 * T ** 2 * X[2:11:4] + T ** 3 * X[3:12:4]/6
    Val[1:10:4] = X[1:10:4] + T * X[2:11:4] + 0.5 * T ** 2 * X[3:12:4]
    Val[2:11:4] = X[2:11:4] + T * X[3:12:4]
    Val[3:12:4] = X[3:12:4]
    Val[12] = X[12] + T * X[13]
    Val[13] = X[13]
    Jacob1 = np.mat([[1, T, 0.5 * T ** 2, T ** 3 / 6], [0, 1, T, 0.5 * T ** 2], [0, 0, 1, T], [0, 0, 0, 1]])
    Jacob2 = np.mat([[1, T], [0, 1]])
    Jacob = block_diag(Jacob1,Jacob1,Jacob1,Jacob2)
    
    return Val, Jacob


def ObservationM(pos, X, Y, Z):
    Val = np.zeros(shape = (len(X), 1))
    dist = np.zeros(shape = (len(X), 1))
    Jacob = np.zeros(shape = (len(X), len(pos)))
    
    for i in range(0, len(X)):
        Val[i] = ((X[i] - pos[0]) ** 2 + (Y[i] - pos[4]) ** 2 + (Z[i] - pos[8]) ** 2) ** 0.5
        dist[i] = Val[i] + pos[12]
        Jacob[i][0] = (pos[0] - X[i])/Val[i]
        Jacob[i][4] = (pos[4] - Y[i])/Val[i]
        Jacob[i][8] = (pos[8] - Z[i])/Val[i]
        Jacob[i][12] = 1
    return dist, Jacob



def efix2topo(X,Y,Z,XYZ_Station):
    
    rtrans = mat([X - XYZ_Station[0], Y - XYZ_Station[1], Z - XYZ_Station[2]])
    rtrans = rtrans.T
    Q1 = mat([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
    
    latitudeXYZ_Station = math.atan2(XYZ_Station[2], math.sqrt(XYZ_Station[0] ** 2 + XYZ_Station[1] ** 2))
    longitudeXYZ_Station = math.atan2(XYZ_Station[1],XYZ_Station[0])
    
    R2 = mat([[math.cos(math.pi / 2 - latitudeXYZ_Station), 0, -math.sin(math.pi / 2 - latitudeXYZ_Station)],\
              [0, 1, 0],\
              [math.sin(math.pi / 2 - latitudeXYZ_Station), 0,  math.cos(math.pi / 2 - latitudeXYZ_Station)]])
    R3 = mat([[math.cos(longitudeXYZ_Station), math.sin(longitudeXYZ_Station), 0],\
              [-math.sin(longitudeXYZ_Station), math.cos(longitudeXYZ_Station),0],\
              [0, 0 ,1]])
    
    r4 = Q1 * R2 * R3 * rtrans
    
    azimuth = math.atan2(float(r4[1]), float(r4[0])) * 180 / math.pi
    elevation = math.atan2(float(r4[2]), math.sqrt(float(r4[0]) ** 2 + float(r4[1]) ** 2)) * 180 / math.pi
    return azimuth, elevation


def Geodetic(r):

    R_equ = 6378.137e3         # Earth's radius [m]; WGS-84
    f     = 1.0/298.257223563  # Flattening; WGS-84   
    eps = 2.2204e-16
    
    epsRequ = eps * R_equ
    e2      = f * (2.0 - f)    # Square of eccentricity
    
    X = r[0]                   # Cartesian coordinates
    Y = r[1]
    Z = r[2]
    rho2 = X ** 2 + Y ** 2     # Square of distance from z-axis
    
#    # Check validity of input data
#    if np.linalg.norm(r) == 0.0:
#       print("invalid input in Geodetic constructor\n")
#       lon = 0.0
#       lat = 0.0
#       h = -R_Earth
#       return
#    end
    
    #Iteration 
    dZ = e2 * Z
    
    while 1:
      ZdZ    =  Z + dZ
      Nh     =  math.sqrt ( rho2 + ZdZ*ZdZ ) 
      SinPhi =  ZdZ / Nh                    # Sine of geodetic latitude
      N      =  R_equ / math.sqrt(1.0-e2*SinPhi*SinPhi)
      dZ_new =  N*e2*SinPhi
      if abs(dZ-dZ_new) < epsRequ:
          break
      dZ = dZ_new
    
    # Longitude, latitude, altitude
    lon = math.atan2 ( Y, X ) * 180 / math.pi
    lat = math.atan2 ( ZdZ, math.sqrt(rho2) ) * 180 / math.pi
    h   = Nh - N
    
    return lon, lat, h


if __name__ == '__main__':

        ephemeris, _, _, _, _, = ReadRinexNav("./meas/Nav/03072019.nav")
        gps, _, _, _ = ReadRinexSprientObs("./meas/Obs/0910circle_point.obs")
        r, NLOSDection = Kalman_Standalone(gps,ephemeris)