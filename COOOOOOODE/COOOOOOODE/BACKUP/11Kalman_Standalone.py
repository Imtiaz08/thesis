



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
    Sf = 36.0; Sg = 1; sigma=5.0;                      
    Qb = np.mat([[Sf * T + Sg * T ** 3 / 3, Sg * T ** 2 / 2], [Sg * T ** 2 / 2, Sg * T]])                            
    Qxyz = sigma ** 2 * np.mat([[T ** 3 / 3, T ** 2 / 2, T, 1], [T ** 2 / 2, T, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]])
    Q = block_diag(Qxyz,Qxyz,Qxyz,Qb)
    
    epochs = [j[1] for i, j in enumerate(obs[0:-1])]
    epoch = list(set(epochs))
    epoch.sort(key=epochs.index)
    r = []
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
        
        
        indx1 = [i for i, j in enumerate(obs[0:-1]) if j[1] == tow]
        # if the car is moving
#        indx2 = [i for i, j in enumerate(obs[0:-1]) if j[1] == tow]
#        if indx1 == [] or indx2 == []:
#            continue
#        pos[1:10:4] = [[float(car[indx2,5])], [float(car[indx2,6])], [float(car[indx2,7])]]                   
#        pos[2:11:4] = [[float(car[indx2,8])], [float(car[indx2,9])], [float(car[indx2,10])]]                  
#        pos[3:12:4] = [[float(car[indx2,11])], [float(car[indx2,12])], [float(car[indx2,13])]]
        # if the car is static
        pos[1:10:4] = [[0], [0], [0]]
        pos[2:11:4] = [[0], [0], [0]]
        pos[3:12:4] = [[0], [0], [0]]
        
        Satclocrr = []; X = []; Y = []; Z = []; PR = []; azim = []; elev = []
        for i in indx1:
            satorbit = ComputeOrbit(obs[i],ephemeris,XYZ_Station)
            if obs[i][6] >= 35:  #Select data with high S/N
                azimuth, elevation = efix2topo(satorbit[4], satorbit[5], satorbit[6], XYZ_Station)
                if elevation >= 10:
                    Satclocrr.append(satorbit[3])
                    X.append(satorbit[4])
                    Y.append(satorbit[5])
                    Z.append(satorbit[6])
                    PR.append(obs[i][3] + satorbit[3])
                    azim.append(azimuth)
                    elev.append(elevation)

        # Set R  variance of measurement error(pseudorange error)
        Rhoerror = 9                                
        R = np.eye(len(X)) * Rhoerror                        
        #positioning using Kalman Filter
        [pos,P] = Extended_KF(Q, R, P, PR, pos, X, Y, Z);
        
        r.append([])
        r[t].append(float(pos[0]))
        r[t].append(float(pos[4]))
        r[t].append(float(pos[8]))
        XYZ_Station = [float(pos[0]), float(pos[4]), float(pos[8])]
    return r
        
        
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
        r = Kalman_Standalone(gps,ephemeris)