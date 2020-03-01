function [GPS_SOW, GPS_Weeks] = Date2GPSTime(year, month, day, hour, minute, second)

gps_week_start = 'January 6 1980 00:00:00';
modnum = 0; % modnum = 0 for no modulo
utcDate = [year, month, day, hour, minute, second];
tmp = mod((datenum(utcDate) - datenum(gps_week_start))/7,modnum); % (Difference in days)/7 = difference in weeks
GPS_Weeks = floor(tmp);
GPS_SOW = round((tmp-GPS_Weeks)*7*24*3600);