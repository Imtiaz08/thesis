

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


def ReadRinexRealObs(filepath):
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
            column = line
            columns.append(column)
    # Extract observation 
    for i, line in enumerate(columns):
        if ">" in line:
            line = line.split(" ")
            while '' in line:
                    line.remove('')
            timestamp = line[1:7]
            GPS_Week, GPS_SOW = Date2GPSTime(timestamp)
        else:
            if line[0] == 'G' or line[0] == 'C' or line[0] == 'E' or line[0] == 'R':
                satellite = []
                satellite.append(GPS_Week)
                satellite.append(GPS_SOW)
                # 1000 for GPS 2000 for Beidou 3000 for Galileo 4000 for Glonass
                if line[0] == 'G':
                    satellite.append(int(line[1:3])+1000)
                elif line[0] == 'C':
                    continue
#                    satellite.append(int(obs[0][1])*10 + int(obs[0][2])+2000)
                elif line[0] == 'E':
                    continue
#                    satellite.append(int(obs[0][1])*10 + int(obs[0][2])+3000)
                elif line[0] == 'R':
                    continue
#                    satellite.append(int(obs[0][1])*10 + int(obs[0][2])+4000)
                if line[5:18].isspace() == True:
                    satellite.append(0)
                else:
                    satellite.append(float(line[5:18]))
                if line[21:33].isspace() == True:
                    satellite.append(0)
                else:
                    satellite.append(float(line[21:33]))
                if line[41:50].isspace() == True:
                    satellite.append(0)
                else:
                    satellite.append(float(line[41:50]))
                if line[59:65].isspace() == True:
                    satellite.append(0)
                else:
                    satellite.append(float(line[59:65]))
                if line[0] == 'G':
                    gps.append(satellite)
                elif line[0] == 'C':
                    beidou.append(satellite)
                elif line[0] == 'E':
                    galileo.append(satellite)
                elif line[0] == 'R':
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
    if round(float(utc[5])) == 60:
        utc[4] = str(int(utc[4]) + 1)
        utc[5] = '0'
    if round(float(utc[4])) == 60:
        utc[3] = str(int(utc[3]) + 1)
        utc[4] = '0'
    date_time = dt(int(utc[0]),int(utc[1]),int(utc[2]),int(utc[3]),int(utc[4]),round(float(utc[5])))
    
    gpsstartdt = 366 + gps_week_start.toordinal() + (gps_week_start - dt.fromordinal(gps_week_start.toordinal())).total_seconds()/(24*60*60)
    utcdt = 366 + date_time.toordinal() + (date_time - dt.fromordinal(date_time.toordinal())).total_seconds()/(24*60*60)
    
    tmp = (utcdt - gpsstartdt)/7
    GPS_Weeks = math.floor(tmp)
    GPS_SOW = round((tmp-GPS_Weeks)*7*24*3600)
    
    return GPS_Weeks, GPS_SOW

