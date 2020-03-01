
function [Obs]=readRinexObsSprient(filePath)
%% This function opens RINEX 3.02 observation files.
% Follows RINEX 3.02 standard. Reads Multiconstellation observables and
% generates an output matrix

%%% ------ Input--- %%%%%%%%%%
%
%   filePath : path to the RINEX observables file
%
%%% ------ Output--- %%%%%%%%%%
%
%   XYZ_station: ECEF coordinates of reference station or point (estimated by receiver)
%
%   observablesHeader: Cell array containing information of observables for
%   each constellation. to look for GPS information type
%   observablesHeader{'G'}
%
%   obs: Matrix containing observables {'week' 'epoch' 'flag' 'prn' 'C1C' 'D1C' 'S1C'}
%       Different for each constellation.
%
% This functions uses functions str2doubleq for quick conversion from
% string to double. It also uses function Date2GPSTime to convert a date to
% TOW.

fprintf('Loading observables...\n\n');
idFile  = fopen(filePath);

generalHeader = {'week', 'epoch', 'flag', 'prn'};

%Initialzie values
measurementsInterval = -1;

%% Read header
while (true)
    line = fgetl(idFile);                                                   % get line
    splitLine = strsplit(line);                                             % Line splited by spaces
    
    if strfind(line,'APPROX POSITION XYZ')                                  % Receiver aprox position
        XYZ_station=real(str2doubleq(strsplit(strtrim(line(1:60)))));
        Obs.XYZ_Station = XYZ_station;
    elseif ~isempty(strfind(line,'SYS / # / OBS TYPES'))                    % Observation types for the different constellations (C1C, D1 and S1 only  )
        constellation = line(1);
        if constellation        == 'G'
            hasGps = 1;
        elseif constellation    == 'R'
            hasGlonass = 1;
        elseif constellation    == 'C'
            hasBeidou = 1;
        end
        
        nObservables = str2doubleq(line(2:7));                               % Number of observables
        observables = splitLine(3:end - 7);                                 % Take the observables only (Not the text regions)
        observables = [generalHeader, observables];                         % Append the header, as the data will be ordered like specified in obsrvables now.
        
        if nObservables >13 %Two line case
            line2 = fgetl(idFile);
            splitLine2 = strsplit(line2);
            observables = [observables, splitLine2(2:end - 7) ];
        end
        
        observablesHeader{uint8(constellation)} = observables;              % Contains all the observables for the constellations.
        %use constellation letter for indexation
        
    elseif strfind(line,'INTERVAL')
        measurementsInterval=str2doubleq(line(5:10));                        % Measurement intervals (Default 1)
        
    elseif strfind(line,'END OF HEADER')
        break;                                                              % End of header loop
    end
end


if measurementsInterval == -1                                               %If interval not set interval = 1
    measurementsInterval = 1;
end
%% Read body

obs = [];                                                                   % Output matrix
epoch = 0;                                                                  % Epoch counter
nObs = 1;

while(~feof(idFile))                                                        % Until end of file
    line = fgetl(idFile);                                                   % Get line
    splitLine = strsplit(line);                                             % Split line by spaces
    
    if strfind(line, '>')                                                   % New epoch
        epoch = epoch + 1;
        %Read time
        year = str2doubleq(splitLine(2));
        month = str2doubleq(splitLine(3));
        day = str2doubleq(splitLine(4));
        hour = str2doubleq(splitLine(5));
        minute = str2doubleq(splitLine(6));
        second = str2doubleq(splitLine(7));
        time = [year, month, day, hour, minute, second];
        time=real(time);
        %[tow,gpsWeek]=Date2GPSTime(time(1),time(2),time(3),time(4)+time(5)/60+time(6)/3600); %Transform date to seconds of week
        [tow,gpsWeek]=Date2GPSTime(time(1),time(2),time(3),time(4),time(5),time(6));
        
        currentEpoch = tow;                                                 % Output
        currentSatellites = str2doubleq(splitLine(9));                      % Satellite information
        currentFlag = str2doubleq(splitLine(8));                            % flag (use/not use)
        
    else
        error 'Loading not correct, satellites skiped'                     % Simple check, it should never jump if using the right rinex version
    end
    
    if currentSatellites == 0
        fprintf('No satellites in epoch');
    end
    
    for i = 1:real(currentSatellites)                                       % Read the epoch satellites
        line = fgetl(idFile);
        constellation = line(1);                                            % First character indicates the constellation
        prn = str2doubleq ([line(2) line(3)]);                               % Satellites PRN number
        
        nObservables = cellfun('length',observablesHeader(uint8(constellation))) - 4; %The header also includes things that are not measurements
        measurementsPosition = (4:16:16*nObservables+4);                    %Vector containing the columns of the measurements. Each 16 columns theres a measurement
        
        if measurementsPosition(end) > length(line)
            measurementsPosition(end) = length(line);       %Correction of a wierd bug
        end
        
        measurementsValue = zeros(1,nObservables); %Initialize vector to store data
        for m = 1:nObservables % Number of observables in the line (Generally 3)
            value = line(measurementsPosition(m):measurementsPosition(m+1)-2); % Column position of measurement. Measurements take 16 columns
            measurementsValue(m) = str2doubleq(value);                      % convert string to double
        end
        
        measurementsValue = real(measurementsValue);                        % output of str2doubleq is imaginary
%         if measurementsValue(1) == 0                                        % if PSR value equals 0
%             continue;                                                       % Skip line (Satellite has no information on L1)
%         end
        switch constellation                                                %Asign constellation based on first char of line
            
            case 'G' %GPS
                prn = prn+1000;
            case 'R'
                prn = prn+2000;
            case 'S'
                prn = prn+3000;
            case 'E'
                prn = prn+4000;
            case 'C'
                prn = prn+5000;
            otherwise
                error 'Unrecognized constellation'                          %Probably 'J' for QZSS
        end
        
        data = [gpsWeek, currentEpoch,currentFlag,prn,measurementsValue,hour,minute,second];   % store data
        obs{nObs} = real(data);
        nObs= nObs+1;
    end
    
    
end

%Convert cell array to matrix. This is necesary to adapt to the rest of the
%alogrithm. Might give problems when different constellations have more
%observables. 

% solution: zeropad the elements in each cell array with a length shorter than the maxlength

maxLengthCell=max(cellfun('size',obs,2));  %finding the longest vector in the cell array
for i=1:length(obs)
    for j=cellfun('size',obs(i),2)+1:maxLengthCell
         obs{i}(j)=0;   % zeropad the elements in each cell array with a length shorter than the maxlength
    end
end
obs = cell2mat(obs');

count1 = 1;
count2 = 1;
count3 = 1;
for i = 1:length(obs)
    if obs(i,4) < 2000
        Obs.gps.data(count1,:) = obs(i,:);
        count1 = count1 + 1;
    elseif obs(i,4) > 2000 && obs(i,4) < 3000
        Obs.glonass.data(count2,:) = obs(i,:);
        count2 = count2 + 1;
    elseif obs(i,4) > 5000
        Obs.beidou.data(count3,:) = obs(i,:);
        count3 = count3 + 1;
    end
end
% GPS
Obs.gps.col.week = 1;
Obs.gps.col.epoch = 2;
Obs.gps.col.flag = 3;
Obs.gps.col.prn = 4;
Obs.gps.col.C1C = 5;
Obs.gps.col.L1C = 6;
Obs.gps.col.D1C = 7;
Obs.gps.col.S1C = 8;
Obs.gps.col.C2L = 9;
Obs.gps.col.L2L = 10;
Obs.gps.col.D2L = 11;
Obs.gps.col.S2L = 12;
% Glonass
Obs.glonass.col.week = 1;
Obs.glonass.col.epoch = 2;
Obs.glonass.col.flag = 3;
Obs.glonass.col.prn = 4;
Obs.glonass.col.C1C = 5;
Obs.glonass.col.L1C = 6;
Obs.glonass.col.D1C = 7;
Obs.glonass.col.S1C = 8;
Obs.glonass.col.C2C = 9;
Obs.glonass.col.L2C = 10;
Obs.glonass.col.D2C = 11;
Obs.glonass.col.S2C = 12;
%Beidou
Obs.beidou.col.week = 1;
Obs.beidou.col.epoch = 2;
Obs.beidou.col.flag = 3;
Obs.beidou.col.prn = 4;
Obs.glonass.col.C1I = 5;
Obs.glonass.col.L1I = 6;
Obs.glonass.col.D1I = 7;
Obs.glonass.col.S1I = 8;

fclose(idFile);
fprintf('Observables loaded successfully\n\n');
end

%function [Obstype] = Obs_types(observablesHeader)
