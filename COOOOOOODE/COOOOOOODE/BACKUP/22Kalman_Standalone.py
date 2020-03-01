



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

def Kalman_Standalone(obs, ephemeris, car):

    XYZ_Station = obs[-1]
    
    # 
    # Set Q state transition variance
#    Sf = 0.01; Sg = 0.01; sigma=0.01;                      
    Sf = 36000000; Sg = 0.01; sigma=1;                      
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
            pos[[3, 7, 11]] = [[0], [0], [0]]                 #Initial velocity
            pos[12] = 0                                       #Clock Bias
            pos[13] = 0                                       #Clock Drift
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
        
        Satclocrr = []; X = []; Y = []; Z = []; PR = []; azim = []; elev = []; PRN = []
        Rhoerror = []; CNR = []
        for i in indx1:
            satorbit = ComputeOrbit(obs[i],ephemeris[0:-2],XYZ_Station)
            if obs[i][6] >= 0: # and obs[i][2] != 10070:  #Select data with high S/N
                azimuth, elevation = efix2topo(satorbit[4], satorbit[5], satorbit[6], XYZ_Station)
                if elevation >= 10:
                    Satclocrr.append(satorbit[3])
                    PRN.append(satorbit[2])
                    X.append(satorbit[4])
                    Y.append(satorbit[5])
                    Z.append(satorbit[6])
                    CNR.append(obs[i][6])
                    # atmosphere correction
#                    T = 2.3 / math.cos(elevation * math.pi / 180)
#                    temp_r = [car[indx2[0]][2], car[indx2[0]][3], car[indx2[0]][4]]
#                    I = Ionosphereklobuchar(temp_r,elevation,azimuth,tow,ephemeris[-2],ephemeris[-1])
                    PR.append(obs[i][3] + satorbit[3])
                    azim.append(azimuth)
                    elev.append(elevation)
                    
              ###################################################################      
#                    sigma = 293 * 0.1 * math.sqrt(0.5) * pow(10, -(obs[i][6]-2*abs(obs[i][6]-53.99))/20)
#                    Rhoerror.append(sigma ** 2)
#                    Rhoerror.append(0.001)
        # Set R  variance of measurement error(pseudorange error)
        Rhoerror = 9
        R = np.eye(len(X)) * Rhoerror                     
# =============================================================================
        #positioning using Kalman Filter
        [pos,P,DOP,WDOP] = Extended_KF(Q, R, P, PR, pos, X, Y, Z);
# =============================================================================
        GDOP = math.sqrt(DOP[0,0]+DOP[1,1]+DOP[2,2]+DOP[3,3])

        WGDOP = math.sqrt(WDOP[0,0]+WDOP[1,1]+WDOP[2,2]+WDOP[3,3])
        
        DOPenu = DOPxyz2enu([car[indx2[0]][2], car[indx2[0]][3], car[indx2[0]][4]] , DOP)
        
        EDOP_list.append(math.sqrt(DOPenu[0,0]))
        NDOP_list.append(math.sqrt(DOPenu[1,1]))
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
    
    K = Pp * H.T * (H * Pp * H.T + R).I
    
    Xp = mat(Xp)
    
    Xo = Xp + K * Y
    
    
    Po = (Ii - K * H) * Pp
    
    # Calculate DOP
    G = H.T[~(H.T == 0).all(1).A1].T
    
    DOP = (G.T * G).I
    
    
    W = np.eye(len(X)) * 1/np.diag(R)
    
    WDOP = (G.T * W * G).I

    
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


def Ionosphereklobuchar(r,elev,azimuth,tow,alfa,beta):

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
    labda = math.atan2 ( Y, X )
    fi = math.atan2 ( ZdZ, math.sqrt(rho2) )
    
    c        =  2.99792458e8             # speed of light
    deg2semi =  1/180                  # degees to semisircles
    semi2rad =  math.pi                       # semisircles to radians
    deg2rad  =  math.pi/180                  # degrees to radians
    
    a = azimuth*deg2rad;                  # asimuth in radians
    e = elev*deg2semi;                    # elevation angle in
                                          # semicircles
    
    psi = 0.0137 / (e+0.11) - 0.022      # Earth Centered angle
    
    lat_i = fi * deg2semi + psi * math.cos(a)     # Subionospheric lat
    
    if lat_i > 0.416:
        lat_i = 0.416;
    elif lat_i < -0.416:
        lat_i = -0.416;
    
                                          # Subionospheric long
    long_i = labda * deg2semi + (psi * math.sin(a) / math.cos(lat_i * semi2rad))
                                          # Geomagnetic latitude
    lat_m = lat_i + 0.064 * math.cos((long_i - 1.617) * semi2rad)
    
    t = 4.32e4 * long_i + tow
    t = t - 86400 * math.floor(t/86400)                    # Seconds of day
    if t > 86400:
        t = t - 86400
    if t < 0:
        t = t + 86400
    
    sF = 1 + 16 * (0.53 - e) ** 3             # Slant factor
    
                                          # Period of model
    PER = beta[0] + beta[1] * lat_m + beta[2] * lat_m ** 2 + beta[3] * lat_m ** 3
    
    if PER < 72000:
        PER = 72000
    
    x = 2 * math.pi * (t - 50400) / PER             # Phase of the model
                                          # (Max at 14.00 =
                                          # 50400 sec local time)
    
                                          # Amplitud of the model
    AMP = alfa[0] + alfa[1] * lat_m + alfa[2] * lat_m ** 2 + alfa[3] * lat_m ** 3
    if AMP < 0:
        AMP = 0
    
                                          # Ionospheric corr.
    if abs(x) > 1.57:
        dIon1 = sF * (5 * 10 ** -9)
    else:
        dIon1 = sF * (5 * 10 ** -9 + AMP * (1 - x * x / 2 + x ** 4 / 24))
    
    dIon1 = c * dIon1

    return dIon1

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
    
#    gps, _, _, _ = ReadRinexRealObs(crf + '/Data/obs/190523_1.obs')
#    ephemeris, _, _, _, _, = ReadRinexRealNav(crf+'/Data/Nav/190523_1.nav')
#    car = pandas.read_csv(crf + '/Data/umt/190523_1.umt', header = None); car = ExtractCar(car,2)

    gps, _, _, _ = ReadRinexRealObs("./meas/Obs/190523_1.obs")
    ephemeris, _, _, _, _, = ReadRinexRealNav("./meas/Nav/190523_1.nav")
    car = pandas.read_csv("./meas/Car/190523_1.umt", header = None); car = ExtractCar(car,2)


    r, GDOP, WGDOP, EDOP, NDOP, clock = Kalman_Standalone(gps,ephemeris,car)
    
    est_clk = mat(list(zip(* clock))[0])    
    plt.figure()
    plt.plot(est_clk.T, label = 'Receiver Clock offset')
    plt.legend(loc='lower right')
    plt.xlabel('Time(sec)')
    plt.ylabel('m')    
    
    
    est_car_x = mat(list(zip(* r))[1])
    est_car_y = mat(list(zip(* r))[2])
    est_car_z = mat(list(zip(* r))[3])    

    ego_car_x = mat(list(zip(* car))[2])
    ego_car_y = mat(list(zip(* car))[3])
    ego_car_z = mat(list(zip(* car))[4])    
    
    
    
    diff = [ego_car_x - est_car_x, ego_car_y - est_car_y, ego_car_z - est_car_z]
    
    
    plt.figure()
    plt.plot(diff[0].T, label = 'x')
    plt.plot(diff[1].T, label = 'y')
    plt.plot(diff[2].T, label = 'z')
    plt.legend(loc='lower right')
    plt.xlabel('Time(sec)')
    plt.ylabel('m')
    
    
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