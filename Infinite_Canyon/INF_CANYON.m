clear;
close all;


%% (North-South) Canyon Geometry
H_wall = 10;  % Height of Wall in m
D_wall_r = 5;  % Distance to (East) right Wall in m
D_wall_l = 10;  % Distance to (West) left Wand in m


%% GSS7000 Log-File zur Konstellation laden

%Log_Path = '.\';
%CLT = readtable([Log_Path 'GSS7000_LogFile_SatConst_V1A1_Example.csv']);  % Read Constellation Log Table (CLT)
CLT = readtable('sat_data_V1A1.csv');
maxTRow = height(CLT);                            % Number of Lines in Log-Table

PRN = unique(CLT.Sat_ID)


Tc = 3000;  % Timestamp (in ms) of logged Constallation to be analyzed   --> e.g. 3000 [ms] corresponds to 5min.

%% Extract Azimuth and Elevation and System (e.g.GPS and Galileo) at Time Tc
TIME = 0;
Tc_i = 0;
for TRow = 1:maxTRow
    TIME = CLT.Time_ms(TRow);
    if TIME == Tc
        Tc_i = Tc_i+1;              % matrix row-th
        Sat_Tc_Vect(1,Tc_i) = CLT.Sat_PRN(TRow);    % Read Satellite Number
        switch string(CLT.Sat_type(TRow))
            case '0'  %'GPS'
                Sat_Tc_Vect(2,Tc_i) = 1;            % GSS7000 Satellite Type Number, 1- GPS, 10 - Galileo
            case '9'  %'GALILEO'
                Sat_Tc_Vect(2,Tc_i) = 10;            % GSS7000 Satellite Type Number, 1- GPS, 10 - Galileo
        end

        Sat_Tc_Vect(3,Tc_i) = CLT.Azimuth(TRow)*180/pi;    % Read Azimuth Angle of Sat.
        
        Sat_Tc_Vect(4,Tc_i) = CLT.Elevation(TRow)*180/pi;    % Read Elevation Angle of Sat.
        
        
        
        % Determine if Sat. is LOS and/or 1st Order Reflection (and Calc Offset Path)
 
        Az = abs(Sat_Tc_Vect(3,Tc_i));
        El = Sat_Tc_Vect(4,Tc_i);

        if Sat_Tc_Vect(3,Tc_i) < 0          % Sat. is West of (North-South) Canyon
            min_El_LOS = atan(H_wall/(D_wall_l/sin(Az*pi/180)))*180/pi;   % min. Elevation-Angle for LOS Reception 

            if El > min_El_LOS
                Sat_Tc_Vect(5,Tc_i) = 1;     % LOS Reception
            else
                Sat_Tc_Vect(5,Tc_i) = 0;     % LOS blocked
            end

            min_El_Refl = atan(sin(Az*pi/180)*H_wall/(2*D_wall_r+D_wall_l))*180/pi;   % min. Elevation Angle of 1st Order Reflection Reception
            max_El_Refl = atan(sin(Az*pi/180)*H_wall/D_wall_r)*180/pi;                % max. Elevation Angle of 1st Order Reflection Reception

            if (El > min_El_Refl) && (El < max_El_Refl)
                Sat_Tc_Vect(6,Tc_i) = 1;     % 1st Order Reflection Reception
                Offset_Refl_Path = sin(Az*pi/180)*cos(El*pi/180)*2*D_wall_r;   % Offset Path Length of Reflection
                Sat_Tc_Vect(7,Tc_i) = Offset_Refl_Path;     % 1st Order Reflection Reception
            else
                Sat_Tc_Vect(6,Tc_i) = 0;     % No 1st Order Reflection
                Sat_Tc_Vect(7,Tc_i) = 0;     % 1st Order Reflection Reception
            end
                Sat_Tc_Vect(8,Tc_i) = min_El_LOS;     % min. Elevation-Angle for LOS Reception
                Sat_Tc_Vect(9,Tc_i) = min_El_Refl;    % min. Elevation Angle of 1st Order Reflection Reception
                Sat_Tc_Vect(10,Tc_i) = max_El_Refl;    % max. Elevation Angle of 1st Order Reflection Reception
        else
            min_El_LOS = atan(H_wall/(D_wall_r/sin(Az*pi/180)))*180/pi;   % min. Elevation-Angle for LOS Reception 

            if El > min_El_LOS
                Sat_Tc_Vect(5,Tc_i) = 1;     % LOS Reception
            else
                Sat_Tc_Vect(5,Tc_i) = 0;     % LOS blocked
            end

            min_El_Refl = atan(sin(Az*pi/180)*H_wall/(2*D_wall_l+D_wall_r))*180/pi;   % min. Elevation Angle of 1st Order Reflection Reception
            max_El_Refl = atan(sin(Az*pi/180)*H_wall/D_wall_l)*180/pi;                % max. Elevation Angle of 1st Order Reflection Reception

            if (El > min_El_Refl) && (El < max_El_Refl)
                Sat_Tc_Vect(6,Tc_i) = 1;     % 1st Order Reflection Reception
                Offset_Refl_Path = sin(Az*pi/180)*cos(El*pi/180)*2*D_wall_l;   % Offset Path Length of Reflection
                Sat_Tc_Vect(7,Tc_i) = Offset_Refl_Path;     % 1st Order Reflection Reception
            else
                Sat_Tc_Vect(6,Tc_i) = 0;     % No 1st Order Reflection
                Sat_Tc_Vect(7,Tc_i) = 0;     % 1st Order Reflection Reception
            end
                Sat_Tc_Vect(8,Tc_i) = min_El_LOS;     % min. Elevation-Angle for LOS Reception
                Sat_Tc_Vect(9,Tc_i) = min_El_Refl;    % min. Elevation Angle of 1st Order Reflection Reception
                Sat_Tc_Vect(10,Tc_i) = max_El_Refl;    % max. Elevation Angle of 1st Order Reflection Reception

        end
        Sat_Tc_Vect(11,Tc_i) = TIME/10;    %  Timestamp in s
    end
    
end


disp(Sat_Tc_Vect(1:10,:));
disp(Sat_Tc_Vect(11,1));



%% Define User Action File (exemplarily for one Time-Stamp, here 5min.)

% check four Cases for every Sat.:

% Case 1: only LOS Reception
%        >> NO UAF Command Line

% Case 2: NO Reception, LOS is blocked and NO 1st Order Reflection Path available
%        GPS:      00:00:00,SWITCH_SAT,V1_A1_VT1,17,2,
%        GALILEO:  00:00:00,SWITCH_SAT,V1_A1_VT10,5,2,
%                                                 ^ Number of Sat.
%                                            ^^^^ Sat.-System: VT1 - GPS, VT10 - GALILEO

% Case 3: LOS + 1st Order Reflection Reception
%        GPS:      00:00:00,SWITCH_SAT,V1_A1_VT1,1,3,1,0,1,2,<41>0,-1,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,,<41>1,-1,10.0,4.9,0.5,0.1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,,,0,0
%                                                      ^ 0 - LOS is available                                                                                       ^^^ offset Path Length in m of the reflected Path                                                                                           
%                                                                                                                                                              ^^^^ relative Attenuation of the reflected Path compared to LOS (10 means the reflected Signal is 10 times weaker than LOS Signal                                                                                          
%        GALILEO:  00:00:00,SWITCH_SAT,V1_A1_VT10,...

% Case 4: NO LOS, ONLY 1st Order Reflection Reception
%        GPS:      00:00:00,SWITCH_SAT,V1_A1_VT1,1,3,1,1,1,2,<41>0,-1,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,,<41>1,-1,10.0,4.9,0.5,0.1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,,,0,0
%                                                      ^ 1 - LOS is not available                                                                                       ^^^ offset Path Length in m of the reflected Path                                                                                           
%                                                                                                                                                              ^^^^ relative Attenuation of the reflected Path compared to LOS (10 means the reflected Signal is 10 times weaker than LOS Signal                                                                                          
%        GALILEO:  00:00:00,SWITCH_SAT,V1_A1_VT10,...

%% exemplarily Implementation of User Action File Generation at a single Timestamp

fid = fopen('UAF_TEST.act','wt');

N_sat = size(Sat_Tc_Vect,2);

for SAT_i = 1 : N_sat
    SAT_NUM = num2str(Sat_Tc_Vect(1,SAT_i));
    SAT_SYS = num2str(Sat_Tc_Vect(2,SAT_i));
    SAT_OFF = num2str(Sat_Tc_Vect(7,SAT_i));
    TIME_STAMP = num2str(Sat_Tc_Vect(11,SAT_i));   % here only one Timestamp was analyzed 
    
    if ((Sat_Tc_Vect(5,SAT_i)==1)&&(Sat_Tc_Vect(6,SAT_i)==0))       % Case 1 only LOS Reception
        fprintf(fid, '\n');
    elseif ((Sat_Tc_Vect(5,SAT_i)==0)&&(Sat_Tc_Vect(6,SAT_i)==0))   % Case 2  NO Reception
        fprintf(fid, [TIME_STAMP ',SWITCH_SAT,V1_A1_VT' SAT_SYS ',' SAT_NUM ',2,']);
        fprintf(fid, '\n');        
    elseif ((Sat_Tc_Vect(5,SAT_i)==1)&&(Sat_Tc_Vect(6,SAT_i)==1))   % Case 3  LOS + 1st Order Reflection Reception
        fprintf(fid, [TIME_STAMP ',SWITCH_SAT,V1_A1_VT' SAT_SYS ',' SAT_NUM ',3,1,0,1,2,<41>0,-1,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,,<41>1,-1,10.0,' SAT_OFF ',0.5,0.1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,,,0,0']);
        fprintf(fid, '\n');                       
    elseif ((Sat_Tc_Vect(5,SAT_i)==0)&&(Sat_Tc_Vect(6,SAT_i)==1))   % Case 4  NO LOS, ONLY 1st Order Reflection Reception
        fprintf(fid, [TIME_STAMP ',SWITCH_SAT,V1_A1_VT' SAT_SYS ',' SAT_NUM ',3,1,1,1,2,<41>0,-1,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,,<41>1,-1,10.0,' SAT_OFF ',0.5,0.1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,,,0,0']);
        fprintf(fid, '\n');                       
    end
    
end

fclose(fid);

