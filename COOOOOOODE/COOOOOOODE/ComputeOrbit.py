

import math
import numpy as np

def ComputeOrbit(Obs,ephemeris,XYZ_station):

    ## Read Data && set constantstim
    ## del the GPSA and GPSB
    
    week = Obs[0]
    time = Obs[1]                               #..... time of observation
    # here only for GPS, if multi-constellation is used, should be changed here
    prn  = Obs[2] - 1000                          #.... PRN number
    svid = [i[0] for i in ephemeris]            #....... satellite PRN
    m0   = [i[2] for i in ephemeris]            #....... mean anomaly at reference time
    a    = [i[3]**2 for i in ephemeris]       #....... semimajor axis
    dn   = [i[4] for i in ephemeris]            #....... mean motion difference
    e    = [i[5] for i in ephemeris]            #....... eccentricity
    w    = [i[6] for i in ephemeris]            #....... argument of perigee
    cuc  = [i[7] for i in ephemeris]            #....... cosine term, arg. of latitude
    cus  = [i[8] for i in ephemeris]            #....... sine term, arg. of latitude
    crc  = [i[9] for i in ephemeris]            #....... cosine term, radius
    crs  = [i[10] for i in ephemeris]           #....... sine term, radius
    i0   = [i[11] for i in ephemeris]           #....... inclination at reference time
    idot = [i[12] for i in ephemeris]           #....... rate of inclination angle
    cic  = [i[13] for i in ephemeris]           #....... cosine term, inclination
    cis  = [i[14] for i in ephemeris]           #....... sine term, inclination
    omg0 = [i[15] for i in ephemeris]           #....... LoAN at weekly epoch
    odot = [i[16] for i in ephemeris]           #....... rate of right ascension
    toe  = [i[17] for i in ephemeris]           #....... time of ephemeris
    af0  = [i[18] for i in ephemeris]           #....... clock offset
    af1  = [i[19] for i in ephemeris]           #....... clock offset rate
    af2  = [i[1] for i in ephemeris]            #....... clock offset accelaration
    tgd  = [i[20] for i in ephemeris]           #....... group time delay 
    
    meu=398600500000000.0       #....... earth's universal gravitational [m^3/s^2]
    odote=7.2921151467e-05      #....... earth's rotation rate (rad/sec)
    lightspeed=299792458.0      #....... speed of light (m/s)
    F=- 4.442807633e-10         #....... Constant, [sec/(meter)^(1/2)]
    
    ## compute positions for single time
    indx1 = [i for i, j in enumerate(svid) if j == prn]
    tsat=time    
    val = []
    for i in indx1:
        val.append(abs(tsat - toe[i]))
    indx2 = val.index(min(val))

    if tsat - toe[indx1[indx2]] < 0:
        if indx2 == 1:
            pass
        else:
            indx2 = indx2 - 1
    indx=indx1[indx2]
    ## Compute sate;llite position
    

    n0 = math.sqrt( meu / a[indx] ** 3)
    t = tsat - toe[indx]
    n = n0 + dn[indx]
    m = m0[indx] + n*t
    m_dot = n
    
    E = kepOrb2E(m,e[indx])
    #Compute relativistic correction term
    dtr = F * e[indx] * math.sqrt(a[indx]) * math.sin(E);
    
    # Compute satellite clock correction
    clkCorr= (af2[indx] * (tsat - toe[indx]) + af1[indx]) * (tsat-toe[indx]) + af0[indx]+ dtr - tgd[indx]
    
    t = t - clkCorr
    # Calculate Velocity
    E_dot = m_dot / (1 - e[indx] * math.cos(E))
    

    v = math.atan2(math.sqrt(1 - e[indx] ** 2) * math.sin(E) , math.cos(E) - e[indx])
    
    v_dot = math.sin(E) * E_dot * (1 + e[indx] * math.cos(v)) / (math.sin(v) * (1 - e[indx] * math.cos(E)))
    
    phi = v + w[indx]
    
    phi_dot = v_dot
    
    du = cus[indx] * math.sin(2 * phi) + cuc[indx] * math.cos(2 * phi)
    dr = crs[indx] * math.sin(2 * phi) + crc[indx] * math.cos(2 * phi)
    di = cis[indx] * math.sin(2 * phi) + cic[indx] * math.cos(2 * phi)
    
    du_dot = 2 * (cus[indx] * math.cos(2 * phi) - cuc[indx] * math.sin(2 * phi)) * phi_dot
    dr_dot = 2 * (crs[indx] * math.cos(2 * phi) - crc[indx] * math.sin(2 * phi)) * phi_dot
    di_dot = 2 * (cis[indx] * math.cos(2 * phi) - cic[indx] * math.sin(2 * phi)) * phi_dot
    
    u = phi + du
    r = a[indx] * (1 - e[indx] * math.cos(E)) + dr

    i = i0[indx] + di + idot[indx] * t

    u_dot = phi_dot + du_dot
    
    r_dot = a[indx] * e[indx] * math.sin(E) * E_dot + dr_dot
    
    i_dot = idot[indx] + di_dot

    
    xp = r * math.cos(u)

    yp = r * math.sin(u)

    xp_dot = r_dot * math.cos(u) - r * math.sin(u) * u_dot

    
    yp_dot = r_dot * math.sin(u) + r * math.cos(u) * u_dot

    
    omg = omg0[indx] + (odot[indx] - odote) * t - odote * toe[indx]

    omg_dot = odot[indx] - odote

    
    xs = xp * math.cos(omg) - yp * math.cos(i) * math.sin(omg)

    ys = xp * math.sin(omg) + yp * math.cos(i) * math.cos(omg)

    zs = yp * math.sin(i)

    vxs = xp_dot * math.cos(omg) - yp_dot * math.cos(i) * math.sin(omg) + yp * math.sin(i) * i_dot * math.sin(omg) - ys * omg_dot
    
    vys = xp_dot * math.sin(omg) + yp_dot * math.cos(i) * math.cos(omg) - yp * math.sin(i) * i_dot * math.cos(omg) + xs * omg_dot
    
    vzs = yp_dot * math.sin(i) + yp * math.cos(i) * i_dot

    
    # # compute the range
    R = [xs - XYZ_station[0], ys - XYZ_station[1], zs - XYZ_station[2]] 

    R = (R[0] ** 2 + R[1] ** 2 + R[2] ** 2 ) ** 0.5

    
    tau = R / lightspeed

    # earth rotation correction
    phi = -odote * tau
    
    corr = np.dot(np.array([[math.cos(phi),-math.sin(phi)],[math.sin(phi),math.cos(phi)]]), np.array([[xs],[ys]]))

    xs = float(corr[0])

    ys = float(corr[1])

    # light Travel time correction
    xs = xs - vxs * tau
    
    ys = ys - vys * tau
    
    zs = zs - vzs * tau

    
    satorbit = []
    satorbit.append(week)
    satorbit.append(time)
    satorbit.append(prn)
    satorbit.append(clkCorr * 299792458)
    satorbit.append(xs)
    satorbit.append(ys)
    satorbit.append(zs)
    satorbit.append(vxs)
    satorbit.append(vys)
    satorbit.append(vzs)
    
    return satorbit
    
def kepOrb2E(M,e):
# Inputs:  - mean anomaly in radians
#          - eccentricity
# Output: Eccentric anomaly
    if -math.pi < M < 0 or M > math.pi:
        E = M - e
    else:
        E = M + e
        
    check=1
    while check > 1e-09:
        E_new = E + (M - E + e * math.sin(E)) / (1 - e * math.cos(E))
        check = abs(E_new - E)
        E = E_new

    
    return E