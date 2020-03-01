



from ReadRinexNav import ReadRinexNav 
from ReadRinexSprientObs import ReadRinexSprientObs
from ComputeOrbit import ComputeOrbit
import math
import numpy as np
from numpy import mat
from scipy.linalg import block_diag

def Kalman_Collaborative(*args):
    T = 1 # positioning interval
    num_car = int((len(args)-1)/3)
        
    # Set Q state transition variance
    Sf = 10; Sg = 0.01; sigma=5.0;                      
    Qb = np.mat([[Sf * T + Sg * T ** 3 / 3, Sg * T ** 2 / 2], [Sg * T ** 2 / 2, Sg * T]])                            
    Qxyz = sigma ** 2 * np.mat([[T ** 3 / 3, T ** 2 / 2, T, 1], [T ** 2 / 2, T, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]])
    Q = block_diag(Qxyz,Qxyz,Qxyz,Qb)
    P = np.eye(14*(num_car))*10;
    
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
            Q = block_diag(Q,Q)
        epochs.append([j[1] for k, j in enumerate(obs[i][0:-1])])
        epoch.append(list(set(epochs[i])))
        epoch[i].sort(key=epochs[i].index)
        length.append(len(epoch[i]))
    localrange = args[-1]
    r = []
    
    
    
    for t, tow in enumerate(epoch[length.index(min(length))]):
        if t == 0:
            # Set initial values of X and P     
            pos_s = []
            for i in range(0, num_car):
                pos_s.append(initial_position(XYZ_Station[i]))
            pos = mat(pos_s)
            pos = np.reshape(pos,(num_car*14,1))
            r.append([])
            r[t].append(float(pos[0]))
            r[t].append(float(pos[4]))
            r[t].append(float(pos[8]))
            continue


        indx = []
        for i in range(0, num_car):
            indx.append([i for i, j in enumerate(obs[i][0:-1]) if j[1] == tow])
            indx.append([i for i, j in enumerate(car[i]) if j[1] == tow])
        indx.append([i for i, j in enumerate(localrange) if j[1] == tow])
        cnt = 0
        for i in range(0, len(indx)):
            if indx[i] == []:
                cnt += 1
        if cnt != 0:
            r.append([])
            r[t].append(r[t-1][0])
            r[t].append(r[t-1][1])
            r[t].append(r[t-1][2])
            continue
        
        for i in range(0, num_car):
            pos[14*i+1:14*i+10:4] = [[float(car[i][indx[i*2+1],5])],[float(car[i][indx[i*2+1],6])],[float(car[i][indx[i*2+1],7])]]
            pos[14*i+2:14*i+11:4] = [[float(car[i][indx[i*2+1],8])],[float(car[i][indx[i*2+1],9])],[float(car[i][indx[i*2+1],10])]]
            pos[14*i+3:14*i+12:4] = [[float(car[i][indx[i*2+1],11])],[float(car[i][indx[i*2+1],12])],[float(car[i][indx[i*2+1],13])]]

        #######################################################################
        Satclocrr = []; X = []; Y = []; Z = []; PR = []; azim = []; elev = []
        num_sat = [0]; threshold = []
        threshold.append(10); threshold.append(20)
        for i in range(0, num_car):
            for j in range(0, len(indx[i*2])):
                satorbit = ComputeOrbit(obs[i][indx[i*2][j]],eph[i],XYZ_Station[i])
                if obs[i][indx[i*2][j]][6] >= 25:  #Select data with high S/N
                    azimuth, elevation = efix2topo(satorbit[4], satorbit[5], satorbit[6], XYZ_Station[i])
                    if elevation >= threshold[i]:
                        Satclocrr.append(satorbit[3])
                        X.append(satorbit[4])
                        Y.append(satorbit[5])
                        Z.append(satorbit[6])
                        PR.append(obs[i][indx[i*2][j]][3] + satorbit[3])
                        azim.append(azimuth)
                        elev.append(elevation)
            num_sat.append(len(X))
        #######################################################################
        for i in range(0, num_car-1):
        ## Static car case
        ## Add range measurements to PR 
            localrange = math.sqrt((float(car[0][indx[i*2+1],2]) - float(car[i+1][indx[i*2+1],2])) ** 2 +
                                   (float(car[0][indx[i*2+1],3]) - float(car[i+1][indx[i*2+1],3])) ** 2 +
                                   (float(car[0][indx[i*2+1],4]) - float(car[i+1][indx[i*2+1],4])) ** 2)
            PR.append(float(localrange[indx[-1], 2]))
                
                
        # Set R  variance of measurement error(pseudorange error)
        Rhoerrorp = 9                                
        Rp = np.eye(len(X)) * Rhoerrorp         
        Rhoerrorr = 2                                
        Rr = np.eye(num_car-1) * Rhoerrorr
        R = block_diag(Rp,Rr)                         
        
        # positioning using Kalman Filter
        [pos,P] = Extended_KF(num_car, num_sat, Q, R, P, PR, pos, X, Y, Z);
        
        # assign the result to r
        r.append([])
        r[t].append(float(pos[0]))
        r[t].append(float(pos[4]))
        r[t].append(float(pos[8]))
        
        # Update XYZ_Station for computing atmosphere correction
        for i in range(0, num_car):
            XYZ_Station[i] = [float(pos[14*i]), float(pos[14*i+4]), float(pos[14*i+8])]
            
    return r
        
        



def Extended_KF(num_car, num_sat, Q, R, P, PR, pos, X, Y, Z):
    Ii = np.eye(len(pos))
    Ii = mat(Ii)
    Xp = np.zeros(shape = (14*num_car, 1))
    fy = np.zeros(shape = (14*num_car, 14*num_car))
    Gr = np.zeros(shape = (num_car-1, 14*num_car))
    dr = []
    Gr_1 = []
    for i in range(0, num_car):
        Xp[14*i:14*(i+1)], _ = CarVelocity(pos[14*i:14*(i+1)] ,1)
        _, fy[14*i:14*(i+1), 14*i:14*(i+1)] = CarVelocity(Xp[14*i:14*(i+1)] ,1)


                


    pos = np.array(pos)
    fy = mat(fy); P = mat(P)
    
    Pp = fy * P * fy.T + Q
    
    gXp = np.zeros(shape = (len(X), 1))
    H = np.zeros(shape = (len(X), 14*num_car))
    for i in range(0, num_car):
        gXp[num_sat[i]:num_sat[i+1]], H[num_sat[i]:num_sat[i+1], 14*i:14*(i+1)] = \
        ObservationM(Xp[14*i:14*(i+1)], X[num_sat[i]:num_sat[i+1]],Y[num_sat[i]:num_sat[i+1]], Z[num_sat[i]:num_sat[i+1]])
    
    for i in range(0, num_car):
        if i < num_car - 1:
        ## Static car case
        ## Add range measurements to PR here
            dr.append(math.sqrt((float(Xp[0]) - float(Xp[(i+1)*14])) ** 2\
                               +(float(Xp[4]) - float(Xp[4+(i+1)*14])) ** 2\
                               +(float(Xp[8]) - float(Xp[8+(i+1)*14])) ** 2))
            Gr_1.append((float(Xp[0]) - float(Xp[(i+1)*14]))/dr[i])
            Gr_1.append((float(Xp[4]) - float(Xp[4+(i+1)*14]))/dr[i])
            Gr_1.append((float(Xp[8]) - float(Xp[8+(i+1)*14]))/dr[i])
            Gr_1.append(0)
            Gr_1 = mat(Gr_1)
            Gr[i, 0:14:4] = Gr_1
            Gr[i, (i+1)*14:(i+2)*14:4] = -Gr_1
    
    
    
    PR = mat(PR); gXp = mat(gXp); dr = mat(dr)
    gXp = np.vstack((gXp, dr))
    
    Y = PR.T - gXp
    
    Pp = mat(Pp); R = mat(R); H = mat(H); Gr = mat(Gr)
    H = np.vstack((H, Gr))
    
    K = Pp * (H.T * np.linalg.inv(H * Pp * H.T + R))
    
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
        dist[i] = Val[i]
        #dist[i] = Val[i] + pos[12]
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



if __name__ == '__main__':
        ephemeris1, _, _, _, _, = ReadRinexNav("C:/Users/lis4hi/Desktop/PythonCode/Data/Nav/03072019.nav")
        gps1, _, _, _ = ReadRinexSprientObs("C:/Users/lis4hi/Desktop/PythonCode/Data/Obs/circle1.obs")
        gps2, _, _, _ = ReadRinexSprientObs("C:/Users/lis4hi/Desktop/PythonCode/Data/Obs/circle2.obs")
        car1 = np.loadtxt("C:/Users/lis4hi/Desktop/PythonCode/Data/Car/circle1.umt")
        car2 = np.loadtxt("C:/Users/lis4hi/Desktop/PythonCode/Data/Car/circle2.umt")
        localrange = np.loadtxt("C:/Users/lis4hi/Desktop/PythonCode/Data/Car/localrange_circle12.txt")
        r = Kalman_Collaborative(gps1,ephemeris1,car1,gps2,ephemeris1,car2,localrange)