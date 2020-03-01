



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
    Sf = 3600000;
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
    for t, tow in enumerate(epoch):
        if t == 0:
            # Set initial values of X and P     
            pos = np.zeros(shape = (14,1))
            
            pos[[0, 4, 8]] = [[XYZ_Station[0]], [XYZ_Station[1]], [XYZ_Station[2]]]               #Initial position
            pos[[1, 5, 9]] = [[0], [0], [0]]
            pos[[2, 6, 10]] = [[0], [0], [0]]
            pos[[3, 7, 11]] = [[0], [0], [0]]                 
            pos[12] = -1000000                                       #Clock Bias
            pos[13] = 150                                  #Clock Drift
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
        
        indx1 = [i for i, j in enumerate(obs) if j[1] == tow]
#         if the car is moving
        indx2 = [i for i, j in enumerate(car) if j[1] == tow]
        if indx1 == [] or indx2 == []:
            continue
        pos[1:10:4] = [[float(car[indx2[0]][5])], [float(car[indx2[0]][6])], [float(car[indx2[0]][7])]]                   
        pos[2:11:4] = [[float(car[indx2[0]][8])], [float(car[indx2[0]][9])], [float(car[indx2[0]][10])]]                  
        pos[3:12:4] = [[0], [0], [0]]
#        # if the car is static
#        pos[1:10:4] = [[0], [0], [0]]
#        pos[2:11:4] = [[0], [0], [0]]
#        pos[3:12:4] = [[0], [0], [0]]
        
        
#        print(t, tow)       # t is the overlap time of both gps and car
#        if t == 2433:
#            pdb.set_trace()                 
#        b = [i for i, val in enumerate(gps) if val[1] == 304798 ]
#        [gps[i] for i in b]


        
        Satclocrr = []; X = []; Y = []; Z = []; PR = []; azim = []; elev = []; PRN = []
        Rhoerror = []; CNR = []
        for i in indx1:
            satorbit = ComputeOrbit(obs[i],ephemeris[0:-2],XYZ_Station)
            if obs[i][6] >= 20 and obs[i][3] > 0: #and obs[i][2] != 10070:  #Select data with high S/N     
                azimuth, elevation = efix2topo(satorbit[4], satorbit[5], satorbit[6], XYZ_Station)
                if elevation >= 5:
                    Satclocrr.append(satorbit[3])
                    PRN.append(satorbit[2])
                    X.append(satorbit[4])
                    Y.append(satorbit[5])
                    Z.append(satorbit[6])
                    CNR.append(obs[i][6])
                    
                     #atmosphere correction
                    T = 2.3 / math.cos((90 - elevation) * math.pi / 180)
                    temp_r = [car[indx2[0]][2], car[indx2[0]][3], car[indx2[0]][4]]
                    I = Ionosphereklobuchar(temp_r,elevation,azimuth,tow,ephemeris[-2],ephemeris[-1])
                    
                    #pdb.set_trace()                    
                    
                    PR.append(obs[i][3] + satorbit[3] - T - I)
                    azim.append(azimuth)
                    elev.append(elevation)
                    
                    sigma = 293 * 0.1 * math.sqrt(0.5) * pow(10, -(obs[i][6]-2*abs(obs[i][6]-53.99))/20)
                    Rhoerror.append(sigma ** 2)
                    #Rhoerror.append(0.001)
        # Set R  variance of measurement error(pseudorange error)
#        Rhoerror = 9
        R = np.eye(len(X)) * Rhoerror       

#        print(t, tow, len(X))
# =============================================================================
        #positioning using Kalman Filter
        if len(X) > 3:
            [pos,P,DOP,WDOP] = Extended_KF(Q, R, P, PR, pos, X, Y, Z);
# =============================================================================
        
            GDOP = math.sqrt(DOP[0,0]+DOP[1,1]+DOP[2,2]+DOP[3,3])
    
            WGDOP = math.sqrt(WDOP[0,0]+WDOP[1,1]+WDOP[2,2]+WDOP[3,3])
            
            DOPenu = DOPxyz2enu([car[indx2[0]][2], car[indx2[0]][3], car[indx2[0]][4]] , DOP)
            
            EDOP_list.append(math.sqrt(DOPenu[0,0]))
            NDOP_list.append(math.sqrt(DOPenu[1,1]))
        else:
            pos, _ = CarVelocity(pos ,1)
            GDOP = np.nan
            WGDOP = np.nan
            EDOP_list.append(np.nan)
            NDOP_list.append(np.nan)
            
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
      
            
            
    return r, GDOP_list, WGDOP_list, EDOP_list, NDOP_list, clock
        
        
def Extended_KF(Q, R, P, PR, pos, X, Y, Z):
    Ii = np.eye(len(pos))
    Ii = mat(Ii)

    
    Xp, fy = CarVelocity(pos ,1)
        
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
    
    K = Pp * (H.T * np.linalg.inv(H * Pp * H.T + R))
    
    Xp = mat(Xp)
    
    Xo = Xp + K * Y
    
    
    Po = (Ii - K * H) * Pp
    
    # Calculate DOP
    G = H.T[~(H.T == 0).all(1).A1].T
    
    DOP =  np.linalg.inv(G.T * G)
    
    
    W = np.eye(len(X)) * 1/np.diag(R)
    
    WDOP =  np.linalg.inv(G.T * W * G)

    return Xo, Po, DOP, WDOP

        
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
        Val[i] = math.sqrt((X[i] - pos[0]) ** 2 + (Y[i] - pos[4]) ** 2 + (Z[i] - pos[8]) ** 2)
        dist[i] = Val[i] + pos[12]
#        dist[i] = Val[i]
        Jacob[i][0] = (pos[0] - X[i])/Val[i]
        Jacob[i][4] = (pos[4] - Y[i])/Val[i]
        Jacob[i][8] = (pos[8] - Z[i])/Val[i]
        Jacob[i][12] = 1
    return dist, Jacob



def Date2GPSTime(utc):
    gps_week_start = dt(1980,1,6,0,0,0)
    date_time = dt(int(utc[0]),int(utc[1]),int(utc[2]),int(utc[3]),int(utc[4]),int(float(utc[5])))
    
    gpsstartdt = 366 + gps_week_start.toordinal() + (gps_week_start - dt.fromordinal(gps_week_start.toordinal())).total_seconds()/(24*60*60)
    utcdt = 366 + date_time.toordinal() + (date_time - dt.fromordinal(date_time.toordinal())).total_seconds()/(24*60*60)
    
    tmp = (utcdt - gpsstartdt)/7
    GPS_Weeks = math.floor(tmp)
    GPS_SOW = round((tmp-GPS_Weeks)*7*24*3600)
    
    return GPS_Weeks, GPS_SOW


def DOPxyz2enu(r,DOP):

    R_equ = 6378.137e3         # Earth's radius [m]; WGS-84
    f     = 1.0/298.257223563  # Flattening; WGS-84   
    eps = 2.2204e-16
    
    epsRequ = eps * R_equ
    e2      = f * (2.0 - f)    # Square of eccentricity
    
    X = r[0]                   # Cartesian coordinates
    Y = r[1]
    Z = r[2]
    rho2 = X ** 2 + Y ** 2     # Square of distance from z-axis

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
    
    # Longitude, latitude, altitude in radian
    lon = math.atan2 ( Y, X )
    lat = math.atan2 ( ZdZ, math.sqrt(rho2) )
    h   = Nh - N
    
    R = np.zeros(shape = (3, 3))
    R[0][0:3] = [-math.sin(lon), -math.cos(lon), 0]
    R[1][0:3] = [-math.cos(lon)*math.sin(lat), -math.sin(lon)*math.sin(lat), math.cos(lat)]
    R[2][0:3] = [math.cos(lon)*math.cos(lat), math.sin(lon)*math.cos(lat), math.sin(lat)]

    DOPenu = R * DOP[0:3, 0:3] * R.T
    
    return DOPenu


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
                #car[i][0] = car[i][0]+6*3600     #########################
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
    #crf = os.getcwd()

#    ephemeris, _, _, _, _, = ReadRinexNav(crf+'/Data/Nav/03072019.nav')
#    gps, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0000.obs')
#    car = pandas.read_csv(crf + '/Data/umt/v2.umt', header = None); car = ExtractCar(car)
    
#    gps, _, _, _ = ReadRinexRealObs("./meas/Obs/190523_1.obs")
#    ephemeris, _, _, _, _, = ReadRinexRealNav("./meas/Nav/190523_1.nav")
#    car = pandas.read_csv("./meas/Car/190523_1.umt", header = None); car = ExtractCar(car,2)
    
    gps, _, _, _ = ReadRinexRealObs("./meas/Obs/newVeh0480.obs")
    ephemeris, _, _, _, _, = ReadRinexRealNav("./meas/Nav/03072019.nav")
    car = pandas.read_csv("./meas/Car/newVeh0480.umt", header = None); car = ExtractCar(car,1)        
    

    r, GDOP, WGDOP, EDOP, NDOP, clock = Kalman_Standalone(gps,ephemeris,car)
    
    
    
    est_clk = np.vstack(list(zip(* clock))[0])    
    clk_dft = np.vstack(list(zip(* clock))[1])
    plt.figure()
    plt.plot(est_clk, label = 'Receiver Clock offset')
    plt.plot(clk_dft, label = 'Receiver Clock drift', color = 'r')
    plt.legend(loc='lower right')
    plt.xlabel('Time(sec)')
    plt.ylabel('m')    
          
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
    print(np.mean(diff_prime[0]), np.mean(diff_prime[1]), np.mean(diff_prime[2]))
    #print(np.mean(diff_prime[0][0, 499:2634]), np.mean(diff_prime[1][0, 1499:2634]), np.mean(diff_prime[2][0, 1499:2634]))

    plt.figure()
    plt.plot()


    t = range(np.size(diff_prime[0]))   
          
    fig, axs = plt.subplots(3, sharex = True)
    fig.suptitle('SPP ENU deviation')
    axs[0].plot(t, diff_prime[0], label = 'xEast')
    axs[1].plot(t, diff_prime[1], label = 'yNorth')
    axs[2].plot(t, diff_prime[2], label = 'zUp')
#    axs[0].set_ylim([-5, 10])
#    axs[1].set_ylim([-10, 20])
#    axs[2].set_ylim([-20, 60])
    axs[0].legend(loc='lower right')
    axs[1].legend(loc='lower right')    
    axs[2].legend(loc='lower right')
    plt.xlabel('Time(sec)')
    plt.ylabel('m')
    
            
    plt.show()
#==============================================================================
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

    
    
    
#    
#    
#    
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