function L = incidentLine(car_pos, sat_azi, sat_ele)
    if (abs(sat_azi)>0 & abs(sat_azi)<pi/2) | sat_azi> pi/2 
        phi = abs(pi/2 - sat_azi);
    else
        phi = 3/2 * pi + sat_azi;
    end

    L = createLine3d(car_pos,pi/2-sat_ele,phi); 
end