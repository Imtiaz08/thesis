function [carmodel] = CarModel(filename,obs)
car = importdata(filename);

if length(car) == 3
    [car(1),car(2),car(3)] = geodetic2ecef(wgs84Ellipsoid('meters'),car(1,1),car(1,2),car(1,3));
    for i = 1:length(obs.gps.data)
        epochsat(i) = obs.gps.data(i,2);
    end
    epochsat = unique(epochsat);
    for i = 1:length(epochsat)
        carmodel(i,1) = obs.gps.data(1,1);
        carmodel(i,2) = epochsat(i);
        carmodel(i,3:5) = car(1:3);
        carmodel(i,6:26) = 0;
    end
else
    for i = 1:length(obs.gps.data)
    epochsat(i) = obs.gps.data(i,13)*3600+obs.gps.data(i,14)*60+obs.gps.data(i,15);
    end
    epochsat = unique(epochsat);
    year = 2019;
    month = 7;
    day = 3;
    epoch = str2num(cell2mat(car.textdata(:,1)));
    count = 1;
    for i = 1:length(epochsat)
        index = find(epoch == epochsat(i));
        if isempty(index)
            continue
        else
            hour = round(epoch(index)/3600);
            minute = round((epoch(index)-3600*hour)/60);
            second = epoch(index)-3600*hour-60*minute;
            time = [year, month, day, hour, minute, second];
            time = real(time);
            [tow,gpsWeek]=Date2GPSTime(time(1),time(2),time(3),time(4),time(5),time(6));
            carmodel(count,1) = gpsWeek;
            carmodel(count,2) = tow;
            carmodel(count,3:26) = car.data(index,:);
            count = count + 1;
        end
    end
end
end


