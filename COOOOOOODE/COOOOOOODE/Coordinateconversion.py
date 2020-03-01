# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 18:14:44 2020

@author: liy4hi
"""
import math
import numpy as np
from numpy import mat


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



def ecef_to_enu(x, y, z, lat0, lon0, h0):
    
    a = 6378137
    b = 6356752.3142
    f = (a - b) / a
    e_sq = f * (2-f)
    
    lamb = math.radians(lat0)
    phi = math.radians(lon0)
    s = math.sin(lamb)
    N = a / math.sqrt(1 - e_sq * s * s)

    sin_lambda = math.sin(lamb)
    cos_lambda = math.cos(lamb)
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)

    x0 = (h0 + N) * cos_lambda * cos_phi
    y0 = (h0 + N) * cos_lambda * sin_phi
    z0 = (h0 + (1 - e_sq) * N) * sin_lambda

    xd = x - x0
    yd = y - y0
    zd = z - z0

    xEast = -sin_phi * xd + cos_phi * yd
    yNorth = -cos_phi * sin_lambda * xd - sin_lambda * sin_phi * yd + cos_lambda * zd
    zUp = cos_lambda * cos_phi * xd + cos_lambda * sin_phi * yd + sin_lambda * zd

    return xEast, yNorth, zUp


def RMSE(predictions, targets):
    predictions = np.array(predictions)
    targets = np.array(targets)
    x = np.sqrt(((predictions - targets) ** 2).mean())
    return x
    
    

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
