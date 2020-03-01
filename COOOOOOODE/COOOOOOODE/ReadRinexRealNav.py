

#=============================================================================
#  D E S C R I P T I O N
#-----------------------------------------------------------------------------
#    Filename: ReadRinexNav.py
#        Name: Shuo Li
#      E-Mail: Shuo.li.chn@gmail.com
#=============================================================================
#  F U N C T I O N
#-----------------------------------------------------------------------------
#   ReadRinexNav Reads a mixed RINEX navigation file *.nav and returns the
#   loaded ephemeris for each constellation
#   Reads Keplerian and Cartesian type ephemeris coming from RINEX 3.x
#=============================================================================
#  E X A M P L E
#-----------------------------------------------------------------------------
#        gps, beidou, galileo, glonass, sbas = ReadRinexNav(filepath)
#       Input: ephemeris filepath 
#      Output: Multiconstellation ephemeris
#-----------------------------------------------------------------------------

__author__ = 'Shuo Li'
__version__ = "1.0"
__email__ = 'Shuo.li.chn@gmail.com'
__status__ = 'Development'


from datetime import datetime as dt
import math


def ReadRinexRealNav(filepath):
	#get navigation messages from GPS/GALILEO/BEIDOU GLONASS
    gps = []
    galileo = []
    glonass = []
    beidou = []
    sbas = []
    
    f = open(filepath, "r")
    fp = filepath
    lookup = "END OF HEADER"
    lookup1 = "GPSA"
    lookup2 = "GPSB"
    for num, line in enumerate(f):
        if lookup1 in line:
            line = line.split(" ")
            line = list(filter(None, line))
            alfa = [float(line[1]), float(line[2]), float(line[3]),float(line[4])]
        if lookup2 in line:
            line = line.split(" ")
            line = list(filter(None, line))
            beta = [float(line[1]), float(line[2]), float(line[3]),float(line[4])]
        if lookup in line:
            index = num


    with open(fp) as fp:
        mylines = fp.readlines()[index + 1:]
        # line_iter = iter(mylines)

        for i in range(0, len(mylines)):
            line = mylines[i]

            if line[0] == "G" or line[0] == "C" or line[0] == "E":
                satellite = []
 		# get navigation msg from GPS/GALILEO/BEIDOU
                line2 = mylines[i + 1]
                line3 = mylines[i + 2]
                line4 = mylines[i + 3]
                line5 = mylines[i + 4]
                line6 = mylines[i + 5]
                line7 = mylines[i + 6]
#                line8 = mylines[i + 7]

                prn = int(line[1:3])
                time = line[4:23]
                af0 = float(line[23:42])
                af1 = float(line[42:61])  # print(line)
                af2 = float(line[61:-1])
                # line_iter.__next__()msg_rinex
#                IODE = float(line2[5:23])
                crs_carrierphase = float(line2[23:42])
                deltan_dopplershift = float(line2[42:61])
                M0 = float(line2[61:-1])
                # line_iter.__next__()
                cuc = float(line3[4:23])
                ecc = float(line3[23:42])
                cus = float(line3[42:61])
                roota = float(line3[61:-1])
                # line_iter.__next__()
                toe = float(line4[4:23])
                cic = float(line4[23:42])
                omega0 = float(line4[42:61])
                cis = float(line4[61:-1])
                # line_iter.__next__()
                i0 = float(line5[4:23])
                crc = float(line5[23:42])
                omega = float(line5[42:61])
                omegadot = float(line5[61:-1])
                # line_iter.__next__()
                idot = float(line6[4:23])
#                codesonl2 = float(line6[23:42])
                week = float(line6[42:61])
                #flag = float(line6[61:-1])
                # line_iter.__next__()
#                svaccuracy = float(line7[4:23])
#                svhealth = float(line7[23:42])
                tgd = float(line7[42:61])
#                iodc = float(line7[61:-1])
                # line_iter.__next__()
#                transtime = float(line8[5:23])
#                fit_interval = float(line8[23:42])
                
                satellite.append(prn)
                satellite.append(af2)
                satellite.append(M0)
                satellite.append(roota)
                satellite.append(deltan_dopplershift)
                satellite.append(ecc)
                satellite.append(omega)
                satellite.append(cuc)
                satellite.append(cus)
                satellite.append(crc)
                satellite.append(crs_carrierphase)
                satellite.append(i0)
                satellite.append(idot)
                satellite.append(cic)
                satellite.append(cis)
                satellite.append(omega0)
                satellite.append(omegadot)
                satellite.append(toe)
                satellite.append(af0)
                satellite.append(af1)
                satellite.append(tgd)

                if line[0] == "G":
                    gps.append(satellite)
                elif line [0] == "E":
                    galileo.append(satellite)
                elif line [0] == "C":
                    beidou.append(satellite)
                
                    
            elif line[0] == "R" or line[0] == "S":
                satellite = []
 		# get navigation msg from GLONASS
                line2 = mylines[i + 1]
                line3 = mylines[i + 2]
                line4 = mylines[i + 3]

                slot_sv = int(line[1:3])
                time = line[4:23]
                time = time.split(" ")
                while '' in time:
                    time.remove('')
                week, toe = Date2GPSTime(time)
                clkbs = float(line[25:42])
                clkrelative = float(line[44:61])  # print(line)
                msgframetime = float(line[63:-1])

                X = float(line2[6:23])
                Xdot = float(line2[25:42])
                Xacc = float(line2[44:61])
                health = float(line2[63:-1])
                # line_iter.__next__()
                Y = float(line3[6:23])
                Ydot = float(line3[25:42])
                Yacc = float(line3[44:61])
                fre_nb = float(line3[63:-1])

                Z = float(line4[6:23])
                Zdot = float(line4[25:42])
                Zacc = float(line4[44:61])
                age = float(line4[63:-1])
                
                satellite.append(slot_sv)
                satellite.append(toe)
                satellite.append(clkbs)
                satellite.append(clkrelative)
                satellite.append(msgframetime)
                satellite.append(X)
                satellite.append(Xdot)
                satellite.append(Xacc)
                satellite.append(health)
                satellite.append(Y)
                satellite.append(Ydot)
                satellite.append(Yacc)
                satellite.append(fre_nb)
                satellite.append(Z)
                satellite.append(Zdot)
                satellite.append(Zacc)
                satellite.append(age)
                satellite.append(1)
                satellite.append(week)
                if line[0] == "R":
                    glonass.append(satellite)
                elif line [0] == "S":
                    sbas.append(satellite)
            else:
#                print("unrecognized constellation")
                continue
        gps.append(alfa)
        gps.append(beta)
    return gps, beidou, galileo, glonass, sbas


def Date2GPSTime(utc):
    gps_week_start = dt(1980,1,6,0,0,0)
    date_time = dt(int(utc[0]),int(utc[1]),int(utc[2]),int(utc[3]),int(utc[4]),int(float(utc[5])))
    
    gpsstartdt = 366 + gps_week_start.toordinal() + (gps_week_start - dt.fromordinal(gps_week_start.toordinal())).total_seconds()/(24*60*60)
    utcdt = 366 + date_time.toordinal() + (date_time - dt.fromordinal(date_time.toordinal())).total_seconds()/(24*60*60)
    
    tmp = (utcdt - gpsstartdt)/7
    GPS_Weeks = math.floor(tmp)
    GPS_SOW = round((tmp-GPS_Weeks)*7*24*3600)
    
    return GPS_Weeks, GPS_SOW