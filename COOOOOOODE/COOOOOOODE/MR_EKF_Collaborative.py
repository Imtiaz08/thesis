
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
from Coordinateconversion import*
from collections import Counter
#args = [gps1,ephemeris1,car1,gps2,ephemeris1,car2,localrange]

def Kalman_Collaborative(*args):
    T = 1 # positioning interval
    num_car = int((len(args))/3)
        
    # Set Q state transition variance
    Sf = 3600000; Sg = 0.01; sigma=1.0;               # sf&sg for receiver clock bias error; sigma affects the convergence sigma(1.0, 5.0)
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
    #localrange = args[-1]
    r = []
    clock = []
    NLOS = {}
    Zglobal = {}
    Zlocal = {}
    test_err = []
    initial_clk = [1200000, 150, -2000000, 150]
    
    
    for t, tow in enumerate(epoch[length.index(min(length))]):
        if t == 0:
            # Set initial values of X and P     
            pos_s = []
            for i in range(0, num_car):
                pos_s.append(initial_position(XYZ_Station[i], initial_clk[2*i], initial_clk[2*i+1]))     # 2*14*1   initial car coordinate of two cars
            pos = mat(pos_s)
            pos = np.reshape(pos,(num_car*14,1))    # 28*1

            
            r.append([])
            clock.append([])
            r[t].append(tow)
            r[t].append(float(pos[0]))
            r[t].append(float(pos[4]))
            r[t].append(float(pos[8]))
            clock[t].append(float(pos[12]))
            clock[t].append(float(pos[13]))
            continue


        indx = []
        for i in range(0, num_car):
            indx.append([i for i, j in enumerate(obs[i][0:-1]) if j[1] == tow])   # make sure tow are the same 
            indx.append([i for i, j in enumerate(car[i]) if j[1] == tow])      
        #indx.append([i for i, j in enumerate(localrange) if j[1] == tow])        # for one certain time step, indx: obs1, car1,obs2,car2,localrange
        cnt = 0
        for i in range(0, len(indx)):
            if indx[i] == []:
                cnt += 1
        if cnt != 0:
#            r.append([])
#            r[t].append(r[t-1][0])        # to see if there is any unmatched meas. between obs,car,local range, if yes, cnt !=0
#            r[t].append(r[t-1][1])        # then apply the position of previous sec to this sec
#            r[t].append(r[t-1][2])
            #r[t].append(r[t-1][3])
            continue
        
        for i in range(0, num_car):      
            pos[14*i+1:14*i+10:4] = [[float(car[i][indx[i*2+1][0]][5])],[float(car[i][indx[i*2+1][0]][6])],[float(car[i][indx[i*2+1][0]][7])]]   # apply the v(vx,vy,vz)
            pos[14*i+2:14*i+11:4] = [[float(car[i][indx[i*2+1][0]][8])],[float(car[i][indx[i*2+1][0]][9])],[float(car[i][indx[i*2+1][0]][10])]]  # apply the a(ax,ay,az)
            pos[14*i+3:14*i+12:4] = [[0], [0], [0]]           # jerk = 0
#            pos[14*i+3:14*i+12:4] = [[float(car[i][indx[i*2+1][0]][11])],[float(car[i][indx[i*2+1][0]][12])],[float(car[i][indx[i*2+1][0]][13])]] # apply the j(jx,jy,jz)

        #########################################################################
        #### num_sat, satclocrr, X,Y,Z, PR, original_PR, azim, ele,SNR ###    dr, Gr 
        Satclocrr = []; X = []; Y = []; Z = []; PR = []; azim = []; elev = []; SNR = []; original_PR = []; PRN = []
       
        num_sat = [0]; threshold = []
        for i in range(num_car+5):
            threshold.append(10)       # threshold of elevation angle for diff car
        for i in range(0, num_car):
            for j in range(0, len(indx[i*2])):
                satorbit = ComputeOrbit(obs[i][indx[i*2][j]],eph[i][0:-2],XYZ_Station[i])  # light travel time, sagnac(sate.pos) relavitistic effect(sate. clk)  
                if obs[i][indx[i*2][j]][6] >= 20 and obs[i][indx[i*2][j]][3] > 0:  #Select data with high S/N
                    azimuth, elevation = efix2topo(satorbit[4], satorbit[5], satorbit[6], XYZ_Station[i])   ## satellite position & receiver position(of last period) 
                    if elevation >= threshold[i]:
                        Satclocrr.append(satorbit[3])
                        PRN.append(satorbit[2])
                        X.append(satorbit[4])            # X Y Z of satellite position 
                        Y.append(satorbit[5])
                        Z.append(satorbit[6])
                        original_PR.append(obs[i][indx[i*2][j]][3]) 
                        SNR.append(obs[i][indx[i*2][j]][6])      
                        
#==============================================================================
 
                        T = 2.3 / math.cos((90 - elevation) * math.pi / 180)
                        temp_r = [car[i][indx[i*2+1][0]][2], car[i][indx[i*2+1][0]][3], car[i][indx[i*2+1][0]][4]]
                        I = Ionosphereklobuchar(temp_r,elevation,azimuth,tow,eph[i][-2],eph[i][-1])
                        
                        PR.append(obs[i][indx[i*2][j]][3] + satorbit[3] - T - I)    # corrected PR measurement with satellite clk error
 
#==============================================================================
                        
                        #PR.append(obs[i][indx[i*2][j]][3] + satorbit[3])    # corrected PR measurement with satellite clk error
                        azim.append(azimuth)
                        elev.append(elevation)
            num_sat.append(len(X))
        #######################################################################
        ## NLOS detection, change the sequence of X,Y,Z,PR, original_PR, azim, elev 
        
        SNR, PR, X, Y, Z, elev, azim, grp1, PRN = SNRsort(SNR, num_sat, num_car, PR, X, Y, Z, elev, azim, PRN)

        
        car_info = []              # at this time, all the car info within network
        for i in range(0, num_car):
            
            netcar = [float(car[i][indx[i*2+1][0]][2]), float(car[i][indx[i*2+1][0]][3]), float(car[i][indx[i*2+1][0]][4])]
            car_info.append(netcar)
        
        
      
        colab_r = []          # localrange w.r.t the ego vehicle
        for i in range(0, num_car-1):
            ## Static car case
            ## Add range measurements to PR 
            colab_range = math.sqrt((float(car[0][indx[1][0]][2]) - float(car[i+1][indx[(i+1)*2+1][0]][2])) ** 2 +
                           (float(car[0][indx[1][0]][3]) - float(car[i+1][indx[(i+1)*2+1][0]][3])) ** 2 +
                           (float(car[0][indx[1][0]][4]) - float(car[i+1][indx[(i+1)*2+1][0]][4])) ** 2)
            colab_r.append(colab_range)
                
        ########################################################################      
        #math.sqrt((car1[0][2] - car2[0][2])**2+(car1[0][3] - car2[0][3])**2+(car1[0][4] - car2[0][4])**2)
        
# =============================================================================
        obs_time = [int(np.mod(tow, 3600*24)/3600), int(np.mod(np.mod(tow, 3600*24), 3600)/60), np.mod(np.mod(np.mod(tow, 3600*24), 3600),60)]
        obs_time_str = ' '.join(str(e) for e in obs_time)
        
        # positioning using Kalman Filter
        
        if len(PR)+num_car-1 >= 4 * num_car:        
            [pos, P, NLOS_t, Zglobal_t, Zlocal_t, err] = Extended_KF(num_car, num_sat, Q, P, PR, pos, X, Y, Z, grp1, car_info, colab_r, PRN);
        
        else:
            for i in range(0, num_car):
                pos[14*i:14*(i+1)], _ = CarVelocity(pos[14*i:14*(i+1)] ,1)
            
            
            
        NLOS[t] = tow, obs_time_str, NLOS_t
        Zglobal[t] = tow, obs_time_str, Zglobal_t
        Zlocal[t] = tow, obs_time_str, Zlocal_t
        test_err.append(err)
        
        # assign the result to r
        r.append([])
        clock.append([])
        r[t].append(tow)
        r[t].append(float(pos[0]))
        r[t].append(float(pos[4]))
        r[t].append(float(pos[8]))
        clock[t].append(float(pos[12]))
        clock[t].append(float(pos[13]))
        
        # Update XYZ_Station for computing atmosphere correction
        for i in range(0, num_car):
            XYZ_Station[i] = [float(pos[14*i]), float(pos[14*i+4]), float(pos[14*i+8])]     #for two cars
        
        
    #return r, Gr
    return r, clock, NLOS, Zglobal, Zlocal, test_err
        
        



def Extended_KF(num_car, num_sat, Q, P, PR, pos, X, Y, Z, grp1, car_info, colab_r, PRN):
    
        
    PRpLRG = PR
    PRpLRG.append(colab_r)
    PR_comb = np.concatenate(PRpLRG)
    X_comb = np.concatenate(X)
    Y_comb = np.concatenate(Y)
    Z_comb = np.concatenate(Z)
    
    
    Ii = np.eye(len(pos))
    Ii = mat(Ii)
    Xp = np.zeros(shape = (14*num_car, 1))
    fy = np.zeros(shape = (14*num_car, 14*num_car))
    dr = []; Gr_1 = []; Gr = np.zeros(shape = (num_car-1, 14*num_car))
    gXp = np.zeros(shape = (len(X_comb), 1))                   #  12*1  estimated range using satellite and receiver postion      
    H = np.zeros(shape = (len(X_comb), 14*num_car))          # 12*28   12: MEAS
    
    for i in range(0, num_car):
        Xp[14*i:14*(i+1)], fy[14*i:14*(i+1), 14*i:14*(i+1)] = CarVelocity(pos[14*i:14*(i+1)] ,1)   # predicted state vector, Fy state transition matrix


    pos = np.array(pos)
    fy = mat(fy); P = mat(P)    
    Pp = fy * P * fy.T + Q                               # predicted covarince matrix 
    
    
    Rweight = np.zeros(shape = (np.size(X_comb), np.size(X_comb)))
    Zg_t = {}    
    Zl_t = {}    
    nlos_t = {}
    
    for i in range(num_car):         # len(nb_car)
        
        othercar = [*range(num_car)]
        othercar.remove(i)
        othercar_range = []
        othercar_pos = []
        for j in othercar:
                localrange = math.sqrt((float(car_info[i][0]) - float(car_info[j][0])) ** 2 +
                                (float(car_info[i][1]) - float(car_info[j][1])) ** 2 +
                                (float(car_info[i][2]) - float(car_info[j][2])) ** 2)
                othercar_range.append(localrange)
                othercar_pos.append(car_info[j])
                
        
        [Gv, Ddb, kesi] = Gv_solution([Xp[14*i], Xp[14*i+4], Xp[14*i+8], Xp[14*i+12]], X[i],Y[i],Z[i], PR[i], othercar_range, othercar_pos)
        #[Gv, Ddb, kesi] = Gv_solution([Xp[14*i], Xp[14*i+4], Xp[14*i+8], 0], X[i],Y[i],Z[i], PR[i], othercar_range, othercar_pos)
        
        nb_current_sate = len(X[i])
        
        Rweight[num_sat[i]:num_sat[i+1],num_sat[i]:num_sat[i+1]], nlos_t_v, Zg_t_v, Zl_t_v, deleteindx = \
        NLOSdect(Gv, Ddb, kesi, grp1[i], num_car, nb_current_sate, PRN[i])    # R of satellites
                                                                                                                # kesi: measured PR
        Zg_t['car' + str(i)] = Zg_t_v
        Zl_t['car' + str(i)] = Zl_t_v
        nlos_t['car' + str(i)] = nlos_t_v
    
    Rhoerrorr = 1                                
    Rr = np.eye(num_car-1) * Rhoerrorr                      # measurement noise of local range
    R_modified = block_diag(Rweight, Rr)                         
    
    

    for i in range(0, num_car):
        gXp[num_sat[i]:num_sat[i+1]], H[num_sat[i]:num_sat[i+1], 14*i:14*(i+1)] = \
        ObservationM(Xp[14*i:14*(i+1)], X_comb[num_sat[i]:num_sat[i+1]],Y_comb[num_sat[i]:num_sat[i+1]], Z_comb[num_sat[i]:num_sat[i+1]])
    
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
    
    PR_comb = mat(PR_comb); gXp = mat(gXp); dr = mat(dr)
    gXp = np.vstack((gXp, dr))         # vertically stack  observation vector
    
    
    Yk = PR_comb.T - gXp                      # .T transpose
    
    Pp = mat(Pp); R_modified = mat(R_modified); H = mat(H); Gr = mat(Gr)
    H = np.vstack((H, Gr))                                 # final geometry matrix
    
    ##########################################
    
   
    K = Pp * (H.T * np.linalg.inv(H * Pp * H.T + R_modified))    # Kalman gain
    
    Xp = mat(Xp)                                      # predicted position 
    
    Xo = Xp + K * Yk                                  # updated position 
    
    
    Po = (Ii - K * H) * Pp                          # updated covariance matrix. Pp predicted cov matrix
    
   
    return Xo, Po, nlos_t, Zg_t, Zl_t, np.mean(Yk)



def SNRsort(SNR, nb_sat, nb_car, corr_meas, X, Y, Z, ele, azim, PRN):
    thre = 42 
    SNR_sep = []
    corr_meas_sep = []
    X_sep = []
    Y_sep = []
    Z_sep = []
    ele_sep = []
    azim_sep = []
    PRN_sep = []
    grp1 = []
    
    SNR_sorted = []
    corr_meas_sorted = []
    X_sorted = []
    Y_sorted = []
    Z_sorted = []
    ele_sorted = []
    azim_sorted = []
    PRN_sorted = []
    
    for i in range(len(nb_sat)-1):                      # nb_car = len(nb_sat)-1
        SNR_sep.append(SNR[nb_sat[i]:nb_sat[i+1]])
        corr_meas_sep.append(corr_meas[nb_sat[i]:nb_sat[i+1]])
        X_sep.append(X[nb_sat[i]:nb_sat[i+1]])
        Y_sep.append(Y[nb_sat[i]:nb_sat[i+1]])
        Z_sep.append(Z[nb_sat[i]:nb_sat[i+1]])
        ele_sep.append(ele[nb_sat[i]:nb_sat[i+1]])
        azim_sep.append(azim[nb_sat[i]:nb_sat[i+1]])
        PRN_sep.append(PRN[nb_sat[i]:nb_sat[i+1]])        
        
    for i in range(nb_car):
        SNR_phi1 = []
        corr_meas_phi1 = []
        X_phi1 = []
        Y_phi1 = []
        Z_phi1 = []
        ele_phi1 = []
        azim_phi1 = []
        PRN_phi1 = []
        
        SNR_phi2 = []
        corr_meas_phi2 = []
        X_phi2 = []
        Y_phi2 = []
        Z_phi2 = []
        ele_phi2 = []
        azim_phi2 = []                
        PRN_phi2 = []        
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
                PRN_phi1.append(PRN_sep[i][j])
                cnt += 1
            else:
                SNR_phi2.append(SNR_sep[i][j])
                corr_meas_phi2.append(corr_meas_sep[i][j])
                X_phi2.append(X_sep[i][j])
                Y_phi2.append(Y_sep[i][j])
                Z_phi2.append(Z_sep[i][j])
                ele_phi2.append(ele_sep[i][j])
                azim_phi2.append(azim_sep[i][j])                       
                PRN_phi2.append(PRN_sep[i][j])
                
        SNR_sorted.append(np.hstack((SNR_phi1,SNR_phi2)))                     # sorted for every single car
        corr_meas_sorted.append(np.hstack((corr_meas_phi1,corr_meas_phi2)))
        X_sorted.append(np.hstack((X_phi1,X_phi2)))
        Y_sorted.append(np.hstack((Y_phi1,Y_phi2)))
        Z_sorted.append(np.hstack((Z_phi1,Z_phi2)))
        ele_sorted.append(np.hstack((ele_phi1,ele_phi2)))
        azim_sorted.append(np.hstack((azim_phi1,azim_phi2)))
        PRN_sorted.append(np.hstack((PRN_phi1, PRN_phi2)))
        grp1.append(cnt)
    
    return SNR_sorted, corr_meas_sorted, X_sorted, Y_sorted, Z_sorted, ele_sorted, azim_sorted, grp1, PRN_sorted
    
    
    
    
#    
#    SNR_sorted_comb = np.concatenate(SNR_sorted)
#    corr_meas_sorted_comb = np.concatenate(corr_meas_sorted)
#    X_sorted_comb = np.concatenate(X_sorted)
#    Y_sorted_comb = np.concatenate(Y_sorted)
#    Z_sorted_comb = np.concatenate(Z_sorted)
#    ele_sorted_comb = np.concatenate(ele_sorted)
#    azim_sorted_comb = np.concatenate(azim_sorted)
#

   
#    
#
#    return R_modified, Z_dect, Z_sgl, nlos_dect, SNR_sorted_comb, corr_meas_sorted_comb, X_sorted_comb,Y_sorted_comb, Z_sorted_comb,ele_sorted_comb,azim_sorted_comb  
#    
    
    
def NLOSdect(Gv, D, meas, L1, num_car, num_sate, PRN):
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
    
    Rhoerrorp = 9                                           
    Rp = np.eye(num_sate) * Rhoerrorp                          # measurement noise of PR
    Rhoerrorr = 1                                
    Rr = np.eye(num_car-1) * Rhoerrorr                      # measurement noise of local range
    R = block_diag(Rp,Rr)                                    # size: (len(X) + num_car - 1) * (len(X) + num_car - 1)    
       
        
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
    
    if L1 > 0 and Z_zeta > 5:
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

    
    
    
    
    
    
def Gv_solution(pre_pos,X,Y,Z, PR, otherrange, otherpos):
    
    dist = np.zeros(shape = (len(X), 1))
    G_sate = np.zeros(shape = (len(X), 4))          #############################################3
    dist_prime = np.zeros(shape = (len(X), 1))      # true range plus clk bias
    nb_range = len(otherrange)
    
    for i in range(0, len(X)):
        dist[i] = np.sqrt((X[i] - pre_pos[0]) ** 2 + (Y[i] - pre_pos[1]) ** 2 + (Z[i] - pre_pos[2]) ** 2)     # sqrt(x^2+y^2+z^2) denominator
        dist_prime[i] = dist[i] + pre_pos[3]            # clk bias 
        G_sate[i][0] = (pre_pos[0] - X[i])/dist[i]
        G_sate[i][1] = (pre_pos[1] - Y[i])/dist[i]
        G_sate[i][2] = (pre_pos[2] - Z[i])/dist[i]
        G_sate[i][3] = 1                             ####################################
    
    lr = []; Gr_1 = []; G_lrg = np.zeros(shape = (nb_range,4))  ##########################
    
    for i in range(0, nb_range):
    ## Add range measurements to PR here
        ego_obs = PR.tolist()            
        ego_obs = ego_obs + otherrange                       # PR add local range to last raw  
        
        lr.append(math.sqrt((float(pre_pos[0]) - float(otherpos[i][0])) ** 2\
                           +(float(pre_pos[1]) - float(otherpos[i][1])) ** 2\
                           +(float(pre_pos[2]) - float(otherpos[i][2])) ** 2))       # sqrt(x^2+y^2+z^2) denominator
        
        
        Gr_1.append((float(pre_pos[0]) - float(otherpos[i][0]))/lr[i])     # (x1-x2)/dr
        Gr_1.append((float(pre_pos[1]) - float(otherpos[i][1]))/lr[i])             # (y1-y2)/dr
        Gr_1.append((float(pre_pos[2]) - float(otherpos[i][2]))/lr[i])            # (z1-z2)/dr
        Gr_1.append(0)              ############################################################
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
    #dist_prime = np.zeros(shape = (len(X), 1))
    
    for i in range(0, len(X)):
        Val[i] = np.sqrt((X[i] - pos[0]) ** 2 + (Y[i] - pos[4]) ** 2 + (Z[i] - pos[8]) ** 2)   # sqrt(x^2+y^2+z^2) denominator
        dist[i] = Val[i] + pos[12]              # pos [12]  clk bias of last time period
        #dist[i] = Val[i]
        #dist_prime[i] = dist[i] + pos[12]        # clk bias of last time period
        Jacob[i][0] = (pos[0] - X[i])/Val[i]
        Jacob[i][4] = (pos[4] - Y[i])/Val[i]
        Jacob[i][8] = (pos[8] - Z[i])/Val[i]
        Jacob[i][12] = 1
    return dist, Jacob




def initial_position(XYZ_Station, clk, dft):
    pos = []
    for i in range(0, 14):
        if i == 0:
            pos.append(XYZ_Station[0])
        elif i == 4:
            pos.append(XYZ_Station[1])
        elif i == 8:
            pos.append(XYZ_Station[2])
        elif i == 12:
            pos.append(clk)
        elif i == 13:
            pos.append(dft)
        else:
            pos.append(0)
    return pos





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
                GPS_SOW = GPS_SOW + 6*3600   ############################################
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
#        ephemeris1, _, _, _, _, = ReadRinexNav("./meas/Nav/03072019.nav")
#        gps1, _, _, _ = ReadRinexSprientObs("./meas/Obs/20190703_Veh0000.obs")
#        gps2, _, _, _ = ReadRinexSprientObs("./meas/Obs/20190703_Veh0001.obs")
#
#
#        car1 = pandas.read_csv("./meas/Car/v0.umt", header = None); car1 = ExtractCar(car1, 1)
#        car2 = pandas.read_csv("./meas/Car/v1.umt", header = None); car2 = ExtractCar(car2, 1)
##            
#        ephemeris1, _, _, _, _, = ReadRinexRealNav("./meas/Nav/190523_1.nav")
#        gps1, _, _, _ = ReadRinexRealObs("./meas/Obs/190523_1.obs")
#        gps2, _, _, _ = ReadRinexRealObs("./meas/Obs/190523_2.obs")        
#        
#        
#        
#        car1 = pandas.read_csv("./meas/Car/190523_1.umt", header = None); car1 = ExtractCar(car1, 2)
#        car2 = pandas.read_csv("./meas/Car/190523_2.umt", header = None); car2 = ExtractCar(car2, 2)


    gps1, _, _, _ = ReadRinexRealObs("./meas/Obs/newVeh0481.obs")
    ephemeris1, _, _, _, _, = ReadRinexRealNav("./meas/Nav/03072019.nav")
    car1 = pandas.read_csv("./meas/Car/newVeh0481.umt", header = None); car1 = ExtractCar(car1,1)       
    
    gps2, _, _, _ = ReadRinexRealObs("./meas/Obs/newVeh0480.obs")
    car2 = pandas.read_csv("./meas/Car/newVeh0480.umt", header = None); car2 = ExtractCar(car2,1)   



    #localrange = np.loadtxt("./meas/Car/localrange_circle12.txt")
    r, clock, NLOSDection, statsglobal, statslocal, test_err = Kalman_Collaborative(gps1,ephemeris1,car1,gps2,ephemeris1,car2)
    
    print(np.mean(test_err))
    
    gg = list(list(zip(* list(statsglobal.values())))[2])     
    ggg = [list(gg[i].values()) for i in range(len(gg))]
    car1_global = np.vstack(list(zip(* ggg[1:]))[0])
    car2_global = np.vstack(list(zip(* ggg[1:]))[1])
    
    car1_NLOS_cnt = sum(i > 5 for i in car1_global)   
    car2_NLOS_cnt = sum(i > 5 for i in car2_global)   
    
    print(car1_NLOS_cnt, car2_NLOS_cnt)
    
    plt.figure()
    plt.plot(car1_global, label = 'Car1 Global test')
    plt.plot(car2_global, label = 'Car2 Global test')
    plt.legend(loc = 'lower right')
    
    
#    car1_txyz = list(zip(* car1))[1:5]
#    stp = car1_txyz[0].index(car2[0][1])
#    endtp = car1_txyz[0].index(car2[-1][1])
#    endtp += 1
#    
#    car_origin = [car1_txyz[1][stp],car1_txyz[2][stp], car1_txyz[3][stp]]
#    lon0, lat0, h0 = Geodetic(car_origin)
            
    est_clk = np.vstack(list(zip(* clock))[0])    
    clk_dft = np.vstack(list(zip(* clock))[1])
    plt.figure()
    plt.plot(est_clk, label = 'Receiver Clock offset')
    plt.plot(clk_dft, label = 'Receiver Clock drift', color = 'r')
    plt.legend(loc='lower right')
    plt.xlabel('Time(sec)')
    plt.ylabel('m')    
        
    
    
#    ego_car_x = mat(car1_txyz[1][stp:endtp])
#    ego_car_y = mat(car1_txyz[2][stp:endtp])
#    ego_car_z = mat(car1_txyz[3][stp:endtp])
#    
#    est_car_x = mat(list(zip(* r))[1])
#    est_car_y = mat(list(zip(* r))[2])
#    est_car_z = mat(list(zip(* r))[3]) 
 
#    xEast1, yNorth1, zUp1 = ecef_to_enu(ego_car_x, ego_car_y, ego_car_z, lat0, lon0, h0)        
#    xEast2, yNorth2, zUp2 = ecef_to_enu(est_car_x, est_car_y, est_car_z, lat0, lon0, h0)    
    

    time = list(set(list(zip(* car1))[1]) & set(list(zip(* r))[0]))    
    
    startpos = list(list(zip(* car1))[1]).index(time[0])
    car_origin = [car1[startpos][2],car1[startpos][3], car1[startpos][4]]       #### static case
    lon0, lat0, h0 = Geodetic(car_origin) 

    est_car = [[val[1], val[2], val[3]] for i, val in enumerate(r) if val[0] in time ]
    ego_car = [[val[2], val[3], val[4]] for i, val in enumerate(car1) if val[1] in time ]
    
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
    axs[0].plot(diff_prime[0], color = 'r', label = 'xEast')
    axs[1].plot(diff_prime[1], color = 'r', label = 'yNorth')
    axs[2].plot(diff_prime[2], color = 'r', label = 'zUp')
#    axs[0].set_ylim([-5, 10])
#    axs[1].set_ylim([-10, 20])
#    axs[2].set_ylim([-20, 60])
    axs[0].legend(loc='upper right')
    axs[1].legend(loc='upper right')    
    axs[2].legend(loc='upper right')
    plt.xlabel('Time(sec)')
    plt.ylabel('m')
        
    plt.show()
#==============================================================================
#         diff = [est_car_x - ego_car_x, est_car_y - ego_car_y, est_car_z - ego_car_z]
#   
#         plt.figure()
#         plt.plot(diff[0].T, label = 'x')
#         plt.plot(diff[1].T, label = 'y')
#         plt.plot(diff[2].T, label = 'z')
#         plt.legend(loc='lower right')
#         plt.xlabel('Time(sec)')
#         plt.ylabel('m')
#==============================================================================
        
        #est_clk = est_clk/(1e8*3)
        
        #plt.interactive(True)


        #np.savetxt('./matlab/values.csv', r)
        #np.savetxt('NLOS',)
#==============================================================================
#         file = open('r.txt', 'w')
#         pickle.dump(r, file)
#         file.close()
#==============================================================================
