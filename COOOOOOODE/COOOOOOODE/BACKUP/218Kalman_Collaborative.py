
    


from ReadRinexNav import ReadRinexNav 
from ReadRinexSprientObs import ReadRinexSprientObs
from ReadRinexSprientObs import Date2GPSTime
from ComputeOrbit import ComputeOrbit
import math
import numpy as np
from numpy import mat
from scipy.linalg import block_diag
import pandas
import matplotlib.pyplot as plt

from ReadRinexRealObs import*
from ReadRinexRealNav import*



def Kalman_Collaborative(*args):
    T = 1 # positioning interval
    num_car = int((len(args))/3)
        
# =============================================================================
    # Set Q state transition variance
    Sf = 3600000; Sg = 0.01; sigma=5;                      
    Qb = np.mat([[Sf * T + Sg * T ** 3 / 3, Sg * T ** 2 / 2], [Sg * T ** 2 / 2, Sg * T]])                            
    Qxyz = sigma ** 2 * np.mat([[T ** 3 / 3, T ** 2 / 2, T, 1], [T ** 2 / 2, T, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]])
    QQ = block_diag(Qxyz,Qxyz,Qxyz,Qb)
    Q = block_diag(Qxyz,Qxyz,Qxyz,Qb)
# =============================================================================
    # Set P covariance Matirx
    P = np.eye(14*(num_car))*1;
# =============================================================================
    # get attribute from input
    obs = []; eph = []; car = []; XYZ_Station = []; epochs = []; epoch = []; length = []
    for i in range(0, num_car):
        cnt1 = 3*i
        cnt2 = 3*i + 1
        cnt3 = 3*i + 2
        obs.append(args[cnt1])
        eph.append(args[cnt2])
        car.append(args[cnt3])
        XYZ_Station.append(obs[i][-1])
        if i < num_car - 1:
            Q = block_diag(Q,QQ)
        epochs.append([j[1] for k, j in enumerate(obs[i][0:-1])])
        epoch.append(list(set(epochs[i])))
        epoch[i].sort(key=epochs[i].index)
        length.append(len(epoch[i]))
    r = []
    clock = []
    CGDOP_list = []
    CWGDOP_list = []
    CEDOP_list = []
    CNDOP_list = []
# =============================================================================
    # main loop
    # Calculate postion if ego vehcile and its covariance Matrix of every epoch
    for t, tow in enumerate(epoch[length.index(min(length))]):
    # if t=0, set initial value (from obs file) of position
        if t == 0:
            pos_s = []
            for i in range(0, num_car):
                pos_s.append(initial_position(XYZ_Station[i]))
            pos = mat(pos_s)
            pos = np.reshape(pos,(num_car*14,1))
            r.append([])
            clock.append([])
            r[t].append(tow)
            r[t].append(float(pos[0]))
            r[t].append(float(pos[4]))
            r[t].append(float(pos[8]))
            clock[t].append(float(pos[12]))
            clock[t].append(float(pos[13]))
            continue
# =============================================================================
    # Get the data form same time epoch in obs file and unt file. if one of
    # them is not exist, skip this time epoch and continue
        indx = []
        for i in range(0, num_car):
            indx.append([i for i, j in enumerate(obs[i][0:-1]) if j[1] == tow])
            indx.append([i for i, j in enumerate(car[i]) if j[1] == tow])
        cnt = 0
        for i in range(0, len(indx)):
            if indx[i] == []:
                cnt += 1
    # if skipped, the position of this time epoch can be predicted or simplely
    # equal to the last time epoch. Useful when there is a lot missing epochs
        if cnt != 0:
            continue
        #     r.append([])
        #     r[t].append(r[t-1][0])
        #     r[t].append(r[t-1][1])
        #     r[t].append(r[t-1][2])
# =============================================================================
    # Assign the correct velocity, acceleration and jerk 
        for i in range(0, num_car):
             pos[14*i+1:14*i+10:4] = [[float(car[i][indx[i*2+1][0]][5])],[float(car[i][indx[i*2+1][0]][6])],[float(car[i][indx[i*2+1][0]][7])]]
             pos[14*i+2:14*i+11:4] = [[float(car[i][indx[i*2+1][0]][8])],[float(car[i][indx[i*2+1][0]][9])],[float(car[i][indx[i*2+1][0]][10])]]
             pos[14*i+3:14*i+12:4] = [[0.0],[0.0],[0.0]]
#             pos[14*i+1:14*i+10:4] = [[0.0],[0.0],[0.0]]
#             pos[14*i+2:14*i+11:4] = [[0.0],[0.0],[0.0]]
#             pos[14*i+3:14*i+12:4] = [[0.0],[0.0],[0.0]]
# =============================================================================
        # Calculate satellies' position X Y Z, get corrected pesudorange PR 
        # calculate satellite clock correction
        Satclocrr = []; X = []; Y = []; Z = []; PR = []; azim = []; elev = []
        dr = []; Rhoerrorp = []; Rhoerrorr = []
        num_sat = [0]; 
        threshold = []
#        threshold.append(1)
#        threshold.append(1)
#        threshold.append(15)
#        threshold.append(15)
#        threshold.append(15)
#        threshold.append(15)
#        threshold.append(15)
#        threshold.append(15)
#        threshold.append(15)
#        threshold.append(15)
#        threshold.append(15)
#        threshold.append(15)
#        threshold.append(15)
        threshold.append(10)
        threshold.append(10)
        threshold.append(10)
#        threshold.append(0)
#        threshold.append(0)
#        threshold.append(0)
#        threshold.append(0)
#        threshold.append(0)
#        threshold.append(0)
#        threshold.append(0)
#        threshold.append(0)
#        threshold.append(0)
#        threshold.append(0)
        for i in range(0, num_car):
            for j in range(0, len(indx[i*2])):
                satorbit = ComputeOrbit(obs[i][indx[i*2][j]],eph[i][0:-2],XYZ_Station[i])
                if obs[i][indx[i*2][j]][6] >= 0:  #Select data with high S/N
                    azimuth, elevation = efix2topo(satorbit[4], satorbit[5], satorbit[6], XYZ_Station[i])
                    if elevation >= threshold[i]:
                        Satclocrr.append(satorbit[3])
                        X.append(satorbit[4])
                        Y.append(satorbit[5])
                        Z.append(satorbit[6])
# =============================================================================
                        # atmosphere correction
                        if i == 0:
                            T = 2.3 / math.cos(elevation * math.pi / 180)
                            temp_r = [car[0][indx[i*2+1][0]][2], car[0][indx[i*2+1][0]][3], car[0][indx[i*2+1][0]][4]]
                            I = Ionosphereklobuchar(temp_r,elevation,azimuth,tow,eph[i][-2],eph[i][-1])
# =============================================================================
                           # PR.append(obs[i][indx[i*2][j]][3] + satorbit[3] - T -I)     ############################### T,I
                            PR.append(obs[i][indx[i*2][j]][3] + satorbit[3])
                        else:
                            PR.append(obs[i][indx[i*2][j]][3] + satorbit[3])
                        azim.append(azimuth)
                        elev.append(elevation)
#                        sigma = 293 * 0.1 * math.sqrt(0.5) * pow(10, -(obs[i][indx[i*2][j]][6]-2*abs(obs[i][indx[i*2][j]][6]-49))/20)
#                        Rhoerrorp.append(sigma ** 2)
#                        Rhoerrorp.append(0.01)
#                        if i == 0:
#                            Rhoerrorp.append(100)
#                        else:
#                            Rhoerrorp.append(293 * 0.1 * math.sqrt(0.5) * pow(10, -(obs[i][indx[i*2][j]][6]-2*abs(obs[i][indx[i*2][j]][6]-53.99))/20))

            num_sat.append(len(X))
# =============================================================================
        # Build geometry matrix of localrange
        for i in range(0, num_car - 1):
            ## Static car case
            ## Add range measurements to PR here
            local_range = ((float(car[0][indx[1][0]][2])-float(car[i+1][indx[(i+1)*2+1][0]][2])) ** 2\
                     + (float(car[0][indx[1][0]][3])-float(car[i+1][indx[(i+1)*2+1][0]][3])) ** 2\
                     + (float(car[0][indx[1][0]][4])-float(car[i+1][indx[(i+1)*2+1][0]][4])) ** 2) ** 0.5
#            Rhoerrorr.append(local_range * 0.001)
            PR.append(local_range)
            

# =============================================================================
        # Set R variance of measurement error(pseudorange and localrange error)
        Rhoerrorp = 9
        Rp = np.eye(len(X)) * Rhoerrorp
        Rhoerrorr = 1
        Rr = np.eye(num_car-1) * Rhoerrorr
        R = block_diag(Rp,Rr)
# =============================================================================
        # positioning using Kalman Filter
        [pos,P,CDOP,CWDOP] = Extended_KF(num_car, num_sat, Q, R, P, PR, pos, X, Y, Z, dr);
# =============================================================================
        CGDOP = math.sqrt(CDOP[0,0] + CDOP[1,1] + CDOP[2,2] + CDOP[3,3])
        CWGDOP = math.sqrt(CWDOP[0,0] + CWDOP[1,1] + CWDOP[2,2] + CWDOP[3,3])
        CDOPenu = DOPxyz2enu([car[0][indx[i*2+1][0]][2], car[0][indx[i*2+1][0]][3], car[0][indx[i*2+1][0]][4]] , CDOP)
        CEDOP_list.append(math.sqrt(CDOPenu[0,0]))
        CNDOP_list.append(math.sqrt(CDOPenu[1,1]))
# =============================================================================
        # assign the result to r
        r.append([])
        clock.append([])
        r[t].append(tow)
        r[t].append(float(pos[0]))
        r[t].append(float(pos[4]))
        r[t].append(float(pos[8]))
        clock[t].append(float(pos[12]))
        clock[t].append(float(pos[13]))
        CGDOP_list.append(CGDOP)
        CWGDOP_list.append(CWGDOP)
#        localrange.append(PR[-1])
        # Update XYZ_Station for computing atmosphere correction
        for i in range(0, num_car):
            XYZ_Station[i] = [float(pos[14*i]), float(pos[14*i+4]), float(pos[14*i+8])]
    return r, CGDOP_list, CWGDOP_list, CEDOP_list, CNDOP_list, clock
        
        



def Extended_KF(num_car, num_sat, Q, R, P, PR, pos, X, Y, Z, dr):
    Ii = np.eye(len(pos))
    Ii = mat(Ii)
    Xp = np.zeros(shape = (14*num_car, 1))
    fy = np.zeros(shape = (14*num_car, 14*num_car))
    Gr = np.zeros(shape = (num_car-1, 14*num_car))

    for i in range(0, num_car):
        Xp[14*i:14*(i+1)], fy[14*i:14*(i+1), 14*i:14*(i+1)] = CarVelocity(pos[14*i:14*(i+1)] ,1)


    pos = np.array(pos)
    fy = mat(fy); P = mat(P)
    
    Pp = fy * P * fy.T + Q
    
    gXp = np.zeros(shape = (len(X), 1))
    H = np.zeros(shape = (len(X), 14*num_car))
    for i in range(0, num_car):
        gXp[num_sat[i]:num_sat[i+1]], H[num_sat[i]:num_sat[i+1], 14*i:14*(i+1)] = \
        ObservationM(Xp[14*i:14*(i+1)], X[num_sat[i]:num_sat[i+1]],Y[num_sat[i]:num_sat[i+1]], Z[num_sat[i]:num_sat[i+1]])
        
    for i in range(0, num_car - 1):
        dr.append(math.sqrt((float(Xp[0]) - float(Xp[(i+1)*14])) ** 2\
                       +(float(Xp[4]) - float(Xp[4+(i+1)*14])) ** 2\
                       +(float(Xp[8]) - float(Xp[8+(i+1)*14])) ** 2))
        Gr_1 = []
        Gr_1.append((float(Xp[0]) - float(Xp[(i+1)*14]))/dr[i])
        Gr_1.append((float(Xp[4]) - float(Xp[4+(i+1)*14]))/dr[i])
        Gr_1.append((float(Xp[8]) - float(Xp[8+(i+1)*14]))/dr[i])
        Gr_1.append(0)
        Gr_11 = mat(Gr_1)
        Gr[i, 0:14:4] = Gr_11
        Gr[i, (i+1)*14:(i+2)*14:4] = -Gr_11
    
    PR = mat(PR); gXp = mat(gXp); dr = mat(dr)
    gXp = np.vstack((gXp, dr.T))
    
    Y = PR.T - gXp
    
    Pp = mat(Pp); R = mat(R); H = mat(H); Gr = mat(Gr)
    H = np.vstack((H, Gr))
    
    K = Pp * H.T * (H * Pp * H.T + R).I
    
    Xp = mat(Xp)
    
    Xo = Xp + K * Y
    
    
    Po = (Ii - K * H) * Pp
    # Calculate CDOP
    G = H.T[~(H.T == 0).all(1).A1].T
    
    W = np.eye(len(X)+len(Gr)) * 1/np.diag(R)

    CDOP = (G.T * G).I
    CWDOP = (G.T * W * G).I
    
    return Xo, Po, CDOP, CWDOP



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
        #dist[i] = Val[i]
        Jacob[i][0] = (pos[0] - X[i])/Val[i]
        Jacob[i][4] = (pos[4] - Y[i])/Val[i]
        Jacob[i][8] = (pos[8] - Z[i])/Val[i]
        Jacob[i][12] = 1
    return dist, Jacob



def initial_position(XYZ_Station):
    pos = []
    for i in range(0, 14):
        if i == 0:
            pos.append(XYZ_Station[0])
        elif i == 4:
            pos.append(XYZ_Station[1])
        elif i == 8:
            pos.append(XYZ_Station[2])
        else:
            pos.append(0)
    return pos


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
#    crf = os.getcwd()
#    # Read ephemreis
#    ephemeris, _, _, _, _, = ReadRinexNav(crf+'/Data/Nav/03072019.nav')
#    # Read Observation
#    gps0, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0000.obs')
#    gps1, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0001.obs')
##    gps2, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0002.obs')
##    gps3, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0003.obs')
##    gps4, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0004.obs')
##    gps5, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0005.obs')
##    gps6, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0006.obs')
##    gps7, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0007.obs')
##    gps8, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0008.obs')
##    gps9, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0009.obs')
##    gps10, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0010.obs')
##    gps11, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0011.obs')
##    gps12, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0012.obs')
##    gps13, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0013.obs')
##    gps14, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/20190703_Veh0014.obs')
##    gps0, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/0923circle_point.obs')
##    gps1, _, _, _ = ReadRinexSprientObs(crf + '/Data/obs/0923circle.obs')
#    # Read car 
#    car0 = pandas.read_csv(crf + '/Data/umt/v0.umt', header = None); car0 = ExtractCar(car0)
#    car1 = pandas.read_csv(crf + '/Data/umt/v1.umt', header = None); car1 = ExtractCar(car1)
##    car2 = pandas.read_csv(crf + '/Data/umt/v2.umt', header = None); car2 = ExtractCar(car2)
##    car3 = pandas.read_csv(crf + '/Data/umt/v3.umt', header = None); car3 = ExtractCar(car3)
##    car4 = pandas.read_csv(crf + '/Data/umt/v4.umt', header = None); car4 = ExtractCar(car4)
##    car5 = pandas.read_csv(crf + '/Data/umt/v5.umt', header = None); car5 = ExtractCar(car5)
##    car6 = pandas.read_csv(crf + '/Data/umt/v6.umt', header = None); car6 = ExtractCar(car6)
##    car7 = pandas.read_csv(crf + '/Data/umt/v7.umt', header = None); car7 = ExtractCar(car7)
##    car8 = pandas.read_csv(crf + '/Data/umt/v8.umt', header = None); car8 = ExtractCar(car8)
##    car9 = pandas.read_csv(crf + '/Data/umt/v9.umt', header = None); car9 = ExtractCar(car9)
##    car10 = pandas.read_csv(crf + '/Data/umt/v10.umt', header = None); car10 = ExtractCar(car10)
##    car11 = pandas.read_csv(crf + '/Data/umt/v11.umt', header = None); car11 = ExtractCar(car11)
##    car12 = pandas.read_csv(crf + '/Data/umt/v12.umt', header = None); car12 = ExtractCar(car12)
##    car13 = pandas.read_csv(crf + '/Data/umt/v13.umt', header = None); car13 = ExtractCar(car13)
##    car14 = pandas.read_csv(crf + '/Data/umt/v14.umt', header = None); car14 = ExtractCar(car14)
##    car0 = pandas.read_csv(crf + '/Data/umt/0923circle_point.umt', header = None); car0 = ExtractCar(car0)
##    car1 = pandas.read_csv(crf + '/Data/umt/0923circle.umt', header = None); car1 = ExtractCar(car1)
#    ##    # copos
##    r,CGDOP,CWGDOP = Kalman_Collaborative(gps0,ephemeris,car0,
##                             gps1,ephemeris,car1,gps2,ephemeris,car2,gps3,ephemeris,car3,gps4,ephemeris,car4,gps5,ephemeris,car5,
##                             gps6,ephemeris,car6,gps7,ephemeris,car7,gps8,ephemeris,car8,gps9,ephemeris,car9,gps10,ephemeris,car10)
#    r,CGDOP,CWGDOP,CEDOP,CNDOP = Kalman_Collaborative(gps0,ephemeris,car0,gps1,ephemeris,car1)
##    r,CGDOP,CWGDOP = Kalman_Collaborative(gps5,ephemeris,car5,gps6,ephemeris,car6)
##    r,CGDOP,CWGDOP = Kalman_Collaborative(gps0,ephemeris,car0,gps1,ephemeris,car1,gps2,ephemeris,car2)
##    r,CGDOP,CWGDOP = Kalman_Collaborative(gps0,ephemeris,car0,gps2,ephemeris,car2,gps3,ephemeris,car3)
##    r,CGDOP,CWGDOP = Kalman_Collaborative(gps0,ephemeris,car0,gps13,ephemeris,car13,gps14,ephemeris,car14)
##    r,CGDOP,CWGDOP = Kalman_Collaborative(gps0,ephemeris,car0,gps2,ephemeris,car2,gps3,ephemeris,car3,gps13,ephemeris,car13,gps14,ephemeris,car14)


    ephemeris1, _, _, _, _, = ReadRinexRealNav("./meas/Nav/190523_1.nav")
    gps1, _, _, _ = ReadRinexRealObs("./meas/Obs/190523_1.obs")
    gps2, _, _, _ = ReadRinexRealObs("./meas/Obs/190523_2.obs")        
    
    
    
    car1 = pandas.read_csv("./meas/Car/190523_1.umt", header = None); car1 = ExtractCar(car1, 2)
    car2 = pandas.read_csv("./meas/Car/190523_2.umt", header = None); car2 = ExtractCar(car2, 2)
#        
    #localrange = np.loadtxt("./meas/Car/localrange_circle12.txt")
    r,CGDOP,CWGDOP,CEDOP,CNDOP, clock = Kalman_Collaborative(gps1,ephemeris1,car1,gps2,ephemeris1,car2)

    car1_txyz = list(zip(* car1))[1:5]
    stp = car1_txyz[0].index(car2[0][1])
    endtp = car1_txyz[0].index(car2[-1][1])
    endtp += 1
    
    ego_car_x = mat(car1_txyz[1][stp:endtp])
    ego_car_y = mat(car1_txyz[2][stp:endtp])
    ego_car_z = mat(car1_txyz[3][stp:endtp])
    
    est_car_x = mat(list(zip(* r))[1])
    est_car_y = mat(list(zip(* r))[2])
    est_car_z = mat(list(zip(* r))[3])
    
    est_clk = mat(list(zip(* clock))[0])
    
    
    diff = [ego_car_x - est_car_x, ego_car_y - est_car_y, ego_car_z - est_car_z]
    
    #plt.interactive(True)
    
    
    
    plt.figure()
    plt.plot(diff[0].T, label = 'x')
    plt.plot(diff[1].T, label = 'y')
    plt.plot(diff[2].T, label = 'z')
    plt.legend(loc='lower right')
    plt.xlabel('Time(sec)')
    plt.ylabel('m')
    
    
    plt.figure()
    plt.plot(est_clk.T, label = 'Receiver Clock offset')
    plt.legend(loc='lower right')
    plt.xlabel('Time(sec)')
    plt.ylabel('m')    
    
    
    
    plt.show()

#
#    file = open('copos.txt', 'w')
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
#    file1 = open('CGDOP.txt', 'w')
#    for num, line in enumerate(CGDOP):
#        file1.write(str(CGDOP[num]))
#        file1.write('\n')
#    file1.close()
#
#    file1 = open('CWGDOP.txt', 'w')
#    for num, line in enumerate(CWGDOP):
#        file1.write(str(CWGDOP[num]))
#        file1.write('\n')
#    file1.close()
#    
#    file1 = open('CEDOP.txt', 'w')
#    for num, line in enumerate(CEDOP):
#        file1.write(str(CEDOP[num]))
#        file1.write('\n')
#    file1.close()
#    
#    file1 = open('CNDOP.txt', 'w')
#    for num, line in enumerate(CNDOP):
#        file1.write(str(CNDOP[num]))
#        file1.write('\n')
#    file1.close()