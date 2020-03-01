

#=============================================================================
#  D E S C R I P T I O N
#-----------------------------------------------------------------------------
#    Filename: ReadRinexSprientObs.py
#        Name: Shuo Li
#      E-Mail: Shuo.li.chn@gmail.com
#=============================================================================
#  F U N C T I O N
#-----------------------------------------------------------------------------
#    Reads Multiconstellation observables and generates an output matrix 
#    Only observation file from Spirent simulator !!!
#    The observation file should be RINEX 3.x format
#=============================================================================
#  E X A M P L E
#-----------------------------------------------------------------------------
#       XYZ_Station, gps, beidou, galileo, glonass = ReadRinexSprientObs(File)
#       Input: Observation filepath 
#      Output: Approx Position XYZ 
#              Multiconstellation Observation
#-----------------------------------------------------------------------------

__author__ = 'Shuo Li'
__version__ = "1.0"
__email__ = 'Shuo.li.chn@gmail.com'
__status__ = 'Development'



from datetime import datetime as dt
import math

#PRN = "PRN"
#C1C = "C1C"
#L1C = "L1C"
#D1C = "D1C"
#S1C = "S1C"
#C2C = "C2C"
#L2C = "L2C"
#D2C = "D2C"
#S2C = "S2C"


def ReadRinexSprientObs(filepath):
    columns = []
    gps = []
    galileo = []
    glonass = []
    beidou = []
    
    f = open(filepath, "r")
    fp = filepath
    
    lookup1 = "APPROX POSITION XYZ"  
    lookup2 = "END OF HEADER"
    # Extract APPROX Position XYZ
    for num, line in enumerate(f):
        if lookup1 in line:
            line = line.split(" ")
            line = list(filter(None, line))
            XYZ_Station = [float(line[0]), float(line[1]), float(line[2])] 
        elif lookup2 in line:
            index = num
            break
    # Find END OF HEADER
    with open(fp, "rt") as fp:
        mylines = fp.readlines()[index+1 :]            #read lines after end of header
        for line in mylines:
            column = line.split(" ")
            column = list(filter(None, column))
            columns.append(column)
    # Extract observation 
    for i, line in enumerate(columns):
        if ">" in line:
            timestamp = columns[i][1:7]
            GPS_Week, GPS_SOW = Date2GPSTime(timestamp)
        else:
            obs = columns[i]
            if obs[0][0] == 'G' or obs[0][0] == 'C' or obs[0][0] == 'E' or obs[0][0] == 'R':
                satellite = []
                satellite.append(GPS_Week)
                satellite.append(GPS_SOW)
                # 1000 for GPS 2000 for Beidou 3000 for Galileo 4000 for Glonass
                if obs[0][0] == 'G':
                    satellite.append(int(obs[0][1])*10 + int(obs[0][2])+1000)
                elif obs[0][0] == 'C':
                    satellite.append(int(obs[0][1])*10 + int(obs[0][2])+2000)
                elif obs[0][0] == 'E':
                    satellite.append(int(obs[0][1])*10 + int(obs[0][2])+3000)
                elif obs[0][0] == 'R':
                    satellite.append(int(obs[0][1])*10 + int(obs[0][2])+4000)
                satellite.append(float(obs[1]))
                satellite.append(float(obs[3]))
                satellite.append(float(obs[5]))
                satellite.append(float(obs[7]))
                satellite.append(float(obs[9]))
                satellite.append(float(obs[11]))
                satellite.append(float(obs[13]))
                satellite.append(float(obs[15]))
                if obs[0][0] == 'G':
                    gps.append(satellite)
                elif obs[0][0] == 'C':
                    beidou.append(satellite)
                elif obs[0][0] == 'E':
                    galileo.append(satellite)
                elif obs[0][0] == 'R':
                    glonass.append(satellite)
            else:
                print("unrecognized constellation")
                continue
    # Assign XYZ_Station from Obs file to the end of the list !!!!!!!!!!!!!
    gps.append(XYZ_Station)
    beidou.append(XYZ_Station)
    galileo.append(XYZ_Station)
    glonass.append(XYZ_Station)
    return gps, beidou, galileo, glonass


def Date2GPSTime(utc):
    gps_week_start = dt(1980,1,6,0,0,0)
    date_time = dt(int(utc[0]),int(utc[1]),int(utc[2]),int(utc[3]),int(utc[4]),int(float(utc[5])))
    
    gpsstartdt = 366 + gps_week_start.toordinal() + (gps_week_start - dt.fromordinal(gps_week_start.toordinal())).total_seconds()/(24*60*60)
    utcdt = 366 + date_time.toordinal() + (date_time - dt.fromordinal(date_time.toordinal())).total_seconds()/(24*60*60)
    
    tmp = (utcdt - gpsstartdt)/7
    GPS_Weeks = math.floor(tmp)
    GPS_SOW = round((tmp-GPS_Weeks)*7*24*3600)
    
    return GPS_Weeks, GPS_SOW

