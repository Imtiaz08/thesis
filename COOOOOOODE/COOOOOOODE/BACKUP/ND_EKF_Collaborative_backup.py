
from ReadRinexNav import ReadRinexNav 
from ReadRinexSprientObs import ReadRinexSprientObs
from ComputeOrbit import ComputeOrbit
import math
import numpy as np
from numpy import mat
from scipy.linalg import block_diag


#args = [gps1,ephemeris1,car1,gps2,ephemeris1,car2,localrange]

def Kalman_Collaborative(*args):
    T = 1 # positioning interval
    num_car = int((len(args)-1)/3)
        
    # Set Q state transition variance
    Sf = 36; Sg = 0.01; sigma=5.0;                      # sf&sg for receiver clock bias error; sigma affects the convergence
    Qb = np.mat([[Sf * T + Sg * T ** 3 / 3, Sg * T ** 2 / 2], [Sg * T ** 2 / 2, Sg * T]])                            
    Qxyz = sigma ** 2 * np.mat([[T ** 3 / 3, T ** 2 / 2, T, 1], [T ** 2 / 2, T, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]])
    Q = block_diag(Qxyz,Qxyz,Qxyz,Qb)
    P = np.eye(14*(num_car))*10;               # P_0  
    
    obs = []; eph = []; car = []; XYZ_Station = []; epochs = []; epoch = []; length = []
    for i in range(0, num_car):
        cnt1 = 3*i
        cnt2 = 3*i + 1
        cnt3 = 3*i + 2
        obs.append(args[cnt1])     #obs of 2 vehicles
        eph.append(args[cnt2])     # ephemeris of 2 vehicles
        car.append(args[cnt3])     # 
        XYZ_Station.append(obs[i][-1])    # XYZ_STATION is the last element of obs
        if i < num_car - 1:
            Q = block_diag(Q,Q)           # process noise matrix
        epochs.append([j[1] for k, j in enumerate(obs[i][0:-1])])   
        epoch.append(list(set(epochs[i])))       # remove the repeated 
        epoch[i].sort(key=epochs[i].index)
        length.append(len(epoch[i]))    
    localrange = args[-1]
    r = []
    NLOS = {}
    
    
    for t, tow in enumerate(epoch[length.index(min(length))]):
        if t == 0:
            # Set initial values of X and P     
            pos_s = []
            for i in range(0, num_car):
                pos_s.append(initial_position(XYZ_Station[i]))     # 2*14*1   initial car coordinate of two cars
            pos = mat(pos_s)
            pos = np.reshape(pos,(num_car*14,1))           # 28*1
            r.append([])
            r[t].append(float(pos[0]))
            r[t].append(float(pos[4]))
            r[t].append(float(pos[8]))
            continue


        indx = []
        for i in range(0, num_car):
            indx.append([i for i, j in enumerate(obs[i][0:-1]) if j[1] == tow])   # make sure tow are the same 
            indx.append([i for i, j in enumerate(car[i]) if j[1] == tow])
        indx.append([i for i, j in enumerate(localrange) if j[1] == tow])        # for one certain time step, indx: obs1, car1,obs2,car2,localrange
        cnt = 0
        for i in range(0, len(indx)):
            if indx[i] == []:
                cnt += 1
        if cnt != 0:
            r.append([])
            r[t].append(r[t-1][0])        # to see if there is any unmatched meas. between obs,car,local range, if yes, cnt !=0
            r[t].append(r[t-1][1])        # then apply the position of previous sec to this sec
            r[t].append(r[t-1][2])
            continue
        
        for i in range(0, num_car):      
            pos[14*i+1:14*i+10:4] = [[float(car[i][indx[i*2+1],5])],[float(car[i][indx[i*2+1],6])],[float(car[i][indx[i*2+1],7])]]   # apply the v(vx,vy,vz)
            pos[14*i+2:14*i+11:4] = [[float(car[i][indx[i*2+1],8])],[float(car[i][indx[i*2+1],9])],[float(car[i][indx[i*2+1],10])]]  # apply the a(ax,ay,az)
            pos[14*i+3:14*i+12:4] = [[float(car[i][indx[i*2+1],11])],[float(car[i][indx[i*2+1],12])],[float(car[i][indx[i*2+1],13])]] # apply the j(jx,jy,jz)

        #########################################################################
        #### num_sat, satclocrr, X,Y,Z, PR, original_PR, azim, ele,SNR ###    dr, Gr 
        Satclocrr = []; X = []; Y = []; Z = []; PR = []; azim = []; elev = []; SNR = []; original_PR = []
       
        num_sat = [0]; threshold = []
        threshold.append(10); threshold.append(20)       # threshold of elevation angle 
        for i in range(0, num_car):
            for j in range(0, len(indx[i*2])):
                satorbit = ComputeOrbit(obs[i][indx[i*2][j]],eph[i],XYZ_Station[i])  # light travel time, sagnac(sate.pos) relavitistic effect(sate. clk)  
                if obs[i][indx[i*2][j]][6] >= 0:  #Select data with high S/N
                    azimuth, elevation = efix2topo(satorbit[4], satorbit[5], satorbit[6], XYZ_Station[i])   ## satellite position & receiver position(of last period) 
                    if elevation >= threshold[i]:
                        Satclocrr.append(satorbit[3])
                        X.append(satorbit[4])            # X Y Z of satellite position 
                        Y.append(satorbit[5])
                        Z.append(satorbit[6])
                        original_PR.append(obs[i][indx[i*2][j]][3]) 
                        SNR.append(obs[i][indx[i*2][j]][6])         
                        PR.append(obs[i][indx[i*2][j]][3] + satorbit[3])    # corrected PR measurement with satellite clk error
                        azim.append(azimuth)
                        elev.append(elevation)
            num_sat.append(len(X))
        #######################################################################
        ## NLOS detection, change the sequence of X,Y,Z,PR, original_PR, azim, elev 
        
        car_info = []
        for i in range(0, num_car-1):
            car_info.append([float(car[i+1][indx[i*2+1],2]), float(car[i+1][indx[i*2+1],3]), float(car[i+1][indx[i*2+1],4])])
        
        
        R, Z_dect_coupled, Z_dect_decpl, NLOS_dct, SNR_sort, org_meas_sort, PR, X, Y, Z, ele_sort, azim_sort =\
        SNRdect(SNR, num_sat, num_car, PR, X, Y, Z, elev, azim, pos, indx, localrange, car_info)

        PR = PR.tolist()
        NLOS[t] = NLOS_dct
        
        
        for i in range(0, num_car-1):
            ## Static car case
            ## Add range measurements to PR 
            localrange = math.sqrt((float(car[0][indx[i*2+1],2]) - float(car[i+1][indx[i*2+1],2])) ** 2 +
                           (float(car[0][indx[i*2+1],3]) - float(car[i+1][indx[i*2+1],3])) ** 2 +
                           (float(car[0][indx[i*2+1],4]) - float(car[i+1][indx[i*2+1],4])) ** 2)
            PR.append(float(localrange[indx[-1], 2]))
           
                
        ########################################################################      

        
# =============================================================================
#         ################## Set R  variance of measurement error(pseudorange error)
#         Rhoerrorp = 9                                           
#         Rp = np.eye(len(X)) * Rhoerrorp                          # measurement noise of PR
#         Rhoerrorr = 2                                
#         Rr = np.eye(num_car-1) * Rhoerrorr                      # measurement noise of local range
#         R = block_diag(Rp,Rr)                                    # size: (len(X) + num_car - 1) * (len(X) + num_car - 1)
# =============================================================================
        
        # positioning using Kalman Filter
        [pos,P] = Extended_KF(num_car, num_sat, Q, R, P, PR, pos, X, Y, Z);
        
        # assign the result to r
        r.append([])
        r[t].append(float(pos[0]))      # x,y,z of ego car 
        r[t].append(float(pos[4]))
        r[t].append(float(pos[8]))
        
        # Update XYZ_Station for computing atmosphere correction
        for i in range(0, num_car):
            XYZ_Station[i] = [float(pos[14*i]), float(pos[14*i+4]), float(pos[14*i+8])]     #for two cars
            
    #return r, Gr
    return r, NLOS
        
        



def Extended_KF(num_car, num_sat, Q, R, P, PR, pos, X, Y, Z):
    Ii = np.eye(len(pos))
    Ii = mat(Ii)
    Xp = np.zeros(shape = (14*num_car, 1))
    fy = np.zeros(shape = (14*num_car, 14*num_car))
    dr = []; Gr_1 = []; Gr = np.zeros(shape = (num_car-1, 14*num_car))
    gXp = np.zeros(shape = (len(X), 1))                   #  12*1  estimated range using satellite and receiver postion      
    kesi_D = np.zeros(shape = (len(X),1))    
    H = np.zeros(shape = (len(X), 14*num_car))          # 12*28   12: MEAS
    
    for i in range(0, num_car):
        Xp[14*i:14*(i+1)], _ = CarVelocity(pos[14*i:14*(i+1)] ,1)                    # predicted state vector 
        _, fy[14*i:14*(i+1), 14*i:14*(i+1)] = CarVelocity(Xp[14*i:14*(i+1)] ,1)    # Fy state transition matrix

    pos = np.array(pos)
    fy = mat(fy); P = mat(P)    
    Pp = fy * P * fy.T + Q                               # predicted covarince matrix 
    
    
    

    for i in range(0, num_car):
        gXp[num_sat[i]:num_sat[i+1]], H[num_sat[i]:num_sat[i+1], 14*i:14*(i+1)], kesi_D[num_sat[i]:num_sat[i+1]] = \
        ObservationM(Xp[14*i:14*(i+1)], X[num_sat[i]:num_sat[i+1]],Y[num_sat[i]:num_sat[i+1]], Z[num_sat[i]:num_sat[i+1]])
    
    for i in range(0, num_car-1):       
        ## Static car case
        ## Add range measurements to PR here
        # predicted state vector of this time period 
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
    gXp = np.vstack((gXp, dr))         # vertically stack  observation vector
    
    kesi_D = np.vstack((kesi_D, dr))         #     
    
    Y = PR.T - gXp                      # .T transpose
    
    Pp = mat(Pp); R = mat(R); H = mat(H); Gr = mat(Gr)
    H = np.vstack((H, Gr))                                 # final geometry matrix
    
    ##########################################
    
   
    K = Pp * (H.T * np.linalg.inv(H * Pp * H.T + R))    # Kalman gain
    
    Xp = mat(Xp)                                      # predicted position 
    
    Xo = Xp + K * Y                                  # updated position 
    
    
    Po = (Ii - K * H) * Pp                          # updated covariance matrix. Pp predicted cov matrix
    
   
    return Xo, Po



def SNRdect(SNR, nb_sat, nb_car, corr_meas, X, Y, Z, ele, azim, pos, indx, localrange, car_info):
    thre = 55 
    SNR_sep = []
    corr_meas_sep = []
    X_sep = []
    Y_sep = []
    Z_sep = []
    ele_sep = []
    azim_sep = []
    grp1 = []
    SNR_sorted = []
    corr_meas_sorted = []
    X_sorted = []
    Y_sorted = []
    Z_sorted = []
    ele_sorted = []
    azim_sorted = []

    
    for i in range(len(nb_sat)-1):                      # nb_car = len(nb_sat)-1
        SNR_sep.append(SNR[nb_sat[i]:nb_sat[i+1]])
        corr_meas_sep.append(corr_meas[nb_sat[i]:nb_sat[i+1]])
        X_sep.append(X[nb_sat[i]:nb_sat[i+1]])
        Y_sep.append(Y[nb_sat[i]:nb_sat[i+1]])
        Z_sep.append(Z[nb_sat[i]:nb_sat[i+1]])
        ele_sep.append(ele[nb_sat[i]:nb_sat[i+1]])
        azim_sep.append(azim[nb_sat[i]:nb_sat[i+1]])
        
    for i in range(nb_car):
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

        for j, k in enumerate(SNR_sep[i]):
            if k < thre:
                SNR_phi1.append(SNR_sep[i][j])
                corr_meas_phi1.append(corr_meas_sep[i][j])
                X_phi1.append(X_sep[i][j])
                Y_phi1.append(Y_sep[i][j])
                Z_phi1.append(Z_sep[i][j])
                ele_phi1.append(ele_sep[i][j])
                azim_phi1.append(azim_sep[i][j])
                
                cnt += 1
            else:
                SNR_phi2.append(SNR_sep[i][j])
                corr_meas_phi2.append(corr_meas_sep[i][j])
                X_phi2.append(X_sep[i][j])
                Y_phi2.append(Y_sep[i][j])
                Z_phi2.append(Z_sep[i][j])
                ele_phi2.append(ele_sep[i][j])
                azim_phi2.append(azim_sep[i][j])                       
                
                
        SNR_sorted.append(np.hstack((SNR_phi1,SNR_phi2)))                     # sorted for every single car
        corr_meas_sorted.append(np.hstack((corr_meas_phi1,corr_meas_phi2)))
        X_sorted.append(np.hstack((X_phi1,X_phi2)))
        Y_sorted.append(np.hstack((Y_phi1,Y_phi2)))
        Z_sorted.append(np.hstack((Z_phi1,Z_phi2)))
        ele_sorted.append(np.hstack((ele_phi1,ele_phi2)))
        azim_sorted.append(np.hstack((azim_phi1,azim_phi2)))
        grp1.append(cnt)
    
    
    SNR_sorted_comb = np.concatenate(SNR_sorted)
    corr_meas_sorted_comb = np.concatenate(corr_meas_sorted)
    X_sorted_comb = np.concatenate(X_sorted)
    Y_sorted_comb = np.concatenate(Y_sorted)
    Z_sorted_comb = np.concatenate(Z_sorted)
    ele_sorted_comb = np.concatenate(ele_sorted)
    azim_sorted_comb = np.concatenate(azim_sorted)


    Rweight = np.zeros(shape = (np.size(X_sorted_comb), np.size(X_sorted_comb)))
    Z_dect = []    
    Z_sgl = []    
    nlos_dect = []
    
    for i in range(nb_car):         # len(nb_car)
        
        [Gv, Ddb, kesi] = Gv_solution(pos, X_sorted[i],Y_sorted[i],Z_sorted[i], indx, corr_meas_sorted[i],localrange, nb_car, car_info)
        
        nb_current_sate = np.size(X_sorted[i])
        
        Rweight[nb_sat[i]:nb_sat[i+1],nb_sat[i]:nb_sat[i+1]], nlos, Z, Z1 = NLOSdect(Gv, Ddb, kesi, grp1[i], nb_car, nb_current_sate)    # R of satellites
        
        Z_dect.append(Z)
        Z_sgl.append(Z1)
        nlos_dect.append(nlos)
    
    Rhoerrorr = 1                                
    Rr = np.eye(nb_car-1) * Rhoerrorr                      # measurement noise of local range
    R_modified = block_diag(Rweight, Rr)                         
    

    return R_modified, Z_dect, Z_sgl, nlos_dect, SNR_sorted_comb, corr_meas_sorted_comb, X_sorted_comb,Y_sorted_comb, Z_sorted_comb,ele_sorted_comb,azim_sorted_comb  
    
    
    
def NLOSdect(Gv, D, meas, L1, num_car, num_sate):
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
    Rhoerrorr = 1                                
    Rr = np.eye(num_car-1) * Rhoerrorr                      # measurement noise of local range
    R = block_diag(Rp,Rr)                                    # size: (len(X) + num_car - 1) * (len(X) + num_car - 1)    
       
        
    delta = np.eye(3)
    r_nlos = 10
    
    kesi_V = np.matmul(np.matmul(Gv, delta), Gv.T) + R            # squared matrix
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
                R[i, i] = Z * R[i, i] 
                NLOS_sig[i] = meas_L1
                
            Z_single.append(Z)
                        
    RR = R[:num_sate, :num_sate]
    
    return RR, NLOS_sig, Z_1, Z_single 

    
    
    
    
    
    
def Gv_solution(pos,X,Y,Z,indx,PR_sorted,localrange,nb_car, car_info):
    
    dist = np.zeros(shape = (len(X), 1))
    G_sate = np.zeros(shape = (len(X), 3))
    dist_prime = np.zeros(shape = (len(X), 1))      # true range plus clk bias
        
    
    for i in range(0, len(X)):
        dist[i] = np.sqrt((X[i] - pos[0]) ** 2 + (Y[i] - pos[4]) ** 2 + (Z[i] - pos[8]) ** 2)   # sqrt(x^2+y^2+z^2) denominator
        dist_prime[i] = dist[i] + pos[12]
        G_sate[i][0] = (pos[0] - X[i])/dist[i]
        G_sate[i][1] = (pos[4] - Y[i])/dist[i]
        G_sate[i][2] = (pos[8] - Z[i])/dist[i]
    
    
    lr = []; Gr_1 = []; G_lrg = np.zeros(shape = (nb_car-1,3))
    for i in range(0, nb_car):
        if i < nb_car - 1:
       
        ## Add range measurements to PR here
            ego_obs = PR_sorted.tolist()            
            ego_obs.append(float(localrange[indx[-1], 2]))                       # PR add local range to last raw  
            
            lr.append(math.sqrt((float(pos[0]) - float(pos[(i+1)*14])) ** 2\
                               +(float(pos[4]) - float(pos[4+(i+1)*14])) ** 2\
                               +(float(pos[8]) - float(pos[8+(i+1)*14])) ** 2))     # sqrt(x^2+y^2+z^2) denominator
            Gr_1.append((float(pos[0]) - float(pos[(i+1)*14]))/lr[i])              # (x1-x2)/dr
            Gr_1.append((float(pos[4]) - float(pos[4+(i+1)*14]))/lr[i])             # (y1-y2)/dr
            Gr_1.append((float(pos[8]) - float(pos[8+(i+1)*14]))/lr[i])            # (z1-z2)/dr
            Gr_1 = mat(Gr_1)                                                    # geometry matrix of local range
            G_lrg[i] = Gr_1                                         # 0,4,8,12 
            
            dist_prime = dist_prime.tolist()
            dist_prime.append(lr)
            
    Gv = np.vstack((G_sate,G_lrg))
    
    return Gv, dist_prime, ego_obs







def CarVelocity(X, T):         # predicted position & state transition jacobian 
    Val = np.zeros(shape = (14,1))
    Val[0:9:4] = X[0:9:4] +  T * X[1:10:4] + 0.5 * T ** 2 * X[2:11:4] + T ** 3 * X[3:12:4]/6    # x
    Val[1:10:4] = X[1:10:4] + T * X[2:11:4] + 0.5 * T ** 2 * X[3:12:4]                # v
    Val[2:11:4] = X[2:11:4] + T * X[3:12:4]                                        # a
    Val[3:12:4] = X[3:12:4]                                                       #j
    Val[12] = X[12] + T * X[13]                                             # clk bias of last time stamp + clk drift of last time stamp * time
    Val[13] = X[13]                                                         # clk drift = updated clk drift of last time stamp
    Jacob1 = np.mat([[1, T, 0.5 * T ** 2, T ** 3 / 6], [0, 1, T, 0.5 * T ** 2], [0, 0, 1, T], [0, 0, 0, 1]])
    Jacob2 = np.mat([[1, T], [0, 1]])
    Jacob = block_diag(Jacob1,Jacob1,Jacob1,Jacob2)
    
    return Val, Jacob


def ObservationM(pos, X, Y, Z):             # pseudorange measurement & measurement jacobian  1 SATE.
    Val = np.zeros(shape = (len(X), 1))
    dist = np.zeros(shape = (len(X), 1))
    Jacob = np.zeros(shape = (len(X), len(pos)))
    dist_prime = np.zeros(shape = (len(X), 1))
    
    for i in range(0, len(X)):
        Val[i] = np.sqrt((X[i] - pos[0]) ** 2 + (Y[i] - pos[4]) ** 2 + (Z[i] - pos[8]) ** 2)   # sqrt(x^2+y^2+z^2) denominator
        #dist[i] = Val[i] + pos[12]              # pos [12] ??????????????????????????????????????
        dist[i] = Val[i]
        dist_prime[i] = dist[i] + pos[12]        # clk bias of last time period
        Jacob[i][0] = (pos[0] - X[i])/Val[i]
        Jacob[i][4] = (pos[4] - Y[i])/Val[i]
        Jacob[i][8] = (pos[8] - Z[i])/Val[i]
        Jacob[i][12] = 1
    return dist, Jacob, dist_prime 



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
        ephemeris1, _, _, _, _, = ReadRinexNav("./meas/Nav/03072019.nav")
        gps1, _, _, _ = ReadRinexSprientObs("./meas/Obs/circle1.obs")
        gps2, _, _, _ = ReadRinexSprientObs("./meas/Obs/circle2.obs")
        car1 = np.loadtxt("./meas/Car/circle1.umt")
        car2 = np.loadtxt("./meas/Car/circle2.umt")
        localrange = np.loadtxt("./meas/Car/localrange_circle12.txt")
        r, NLOSDection = Kalman_Collaborative(gps1,ephemeris1,car1,gps2,ephemeris1,car2,localrange)
        np.savetxt('./matlab/values.csv', r)
        #np.savetxt('NLOS',)
#==============================================================================
#         file = open('r.txt', 'w')
#         pickle.dump(r, file)
#         file.close()
#==============================================================================
