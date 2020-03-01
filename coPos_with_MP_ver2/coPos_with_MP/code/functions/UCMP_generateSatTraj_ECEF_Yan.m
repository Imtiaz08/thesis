% SHAIK RIYAAZ UDDIEN- SAPIEZA UNIVERSITY OF ROME, ITALY
% Matlab Code to find out the position of Satellites using Orbital Parameters for time duration of 24 hours
% Satellite Constellation
% Ref: YUMA ALMANAC - US COAST GUARD SITE
function [Xs, Ys, Zs, Xv, Yv, Zv] = UCMP_generateSatTraj_ECEF_Yan(startTime, sampTime, dur, ephemeris, XYZ_station)
meu = 3.986005e14;         % earth's universal gravitational [m^3/s^2]
odote = 7.2921151467e-5;   % earth's rotation rate (rad/sec)
lightspeed = 2.99792458e8; % speed of light (m/s)
F = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]
% Read ephemeris data
m0   = ephemeris.gps(:,3);        %....... mean anomaly at reference time
a    = (ephemeris.gps(:,4)).^2;   %....... semimajor axis
dn   = ephemeris.gps(:,5);        %....... mean motion difference
e    = ephemeris.gps(:,6);        %....... eccentricity
w    = ephemeris.gps(:,7);        %....... argument of perigee
cuc  = ephemeris.gps(:,8);        %....... cosine term, arg. of latitude
cus  = ephemeris.gps(:,9);        %....... sine term, arg. of latitude
crc  = ephemeris.gps(:,10);       %....... cosine term, radius
crs  = ephemeris.gps(:,11);       %....... sine term, radius
i0   = ephemeris.gps(:,12);       %....... inclination at reference time
idot = ephemeris.gps(:,13);       %....... rate of inclination angle
cic  = ephemeris.gps(:,14);       %....... cosine term, inclination
cis  = ephemeris.gps(:,15);       %....... sine term, inclination
omg0 = ephemeris.gps(:,16);       %....... LoAN at weekly epoch
odot = ephemeris.gps(:,17);       %....... rate of right ascension
toe  = ephemeris.gps(:,18);       %....... time of ephemeris
af0  = ephemeris.gps(:,19);       %....... clock offset
af1  = ephemeris.gps(:,20);       %....... clock offset rate
af2  = ephemeris.gps(:,2);        %....... clock offset accelaration
tgd  = ephemeris.gps(:,22);

% Calculate position
numofsat=1:1:32;
for ii=1:length(numofsat)
    % find the certain satellite index
    indx1 = find(ephemeris.gps(:,1) == ii);
    % find the nearest time index
    J= startTime:sampTime:(startTime + dur);
    for ti=1:length(J)
        tsat = J(ti);
        hour = round(tsat/3600);
        minute = round((tsat-3600*hour)/60);
        second = tsat-3600*hour-60*minute;
        [GPS_SOW, ~] = Date2GPSTime(2019, 7, 3,hour, minute, second);
        [~,indx2] = min(abs(ephemeris.gps(indx1(:),18) - GPS_SOW));
        index = indx1(indx2); 
        %%
        n0 = sqrt(meu/a(index)^3);
        t = tsat-toe(index);
        n = n0 + dn(index);
        m = m0(index) + n*t;

        m_dot=n;                                  % Calculate Velocity

        E = kepOrb2E(m,e(index));

        %Compute relativistic correction term
        dtr = F * e(index) * sqrt(a(index)) * sin(E);

        % Compute satellite clock correction    
        clkCorr= (af2(index) * (tsat-toe(index)) + af1(index)) * (tsat-toe(index)) + ...
            af0(index) + dtr - tgd(index);

        t = t - clkCorr;

        E_dot=m_dot/(1-e(index)*cos(E));          % Calculate Velocity    

        v = atan2(sqrt(1-e(index)^2)*sin(E),cos(E)-e(index));

        v_dot=sin(E)*E_dot*(1+e(index)*cos(v))/(sin(v)*(1-e(index)*cos(E)));   % Calculate Velocity

        phi = v + w(index);

        phi_dot=v_dot;                            % Calculate Velocity

        du = cus(index)*sin(2*phi) + cuc(index)*cos(2*phi);
        dr = crs(index)*sin(2*phi) + crc(index)*cos(2*phi);
        di = cis(index)*sin(2*phi) + cic(index)*cos(2*phi);

        du_dot=2*(cus(index)*cos(2*phi)-cuc(index)*sin(2*phi))*phi_dot; % Calculate Velocity
        dr_dot=2*(crs(index)*cos(2*phi)-crc(index)*sin(2*phi))*phi_dot; % Calculate Velocity
        di_dot=2*(cis(index)*cos(2*phi)-cic(index)*sin(2*phi))*phi_dot; % Calculate Velocity

        u = phi + du;
        r = a(index)*(1-e(index)*cos(E)) + dr;
        i = i0(index) + di + idot(index)*t;

        u_dot=phi_dot+du_dot;                         % Calculate Velocity
        r_dot=a(index)*e(index)*sin(E)*E_dot+dr_dot;  % Calculate Velocity
        i_dot=idot(index)+di_dot;                     % Calculate Velocity

        xp = r*cos(u);
        yp = r*sin(u);

        xp_dot=r_dot*cos(u)-r*sin(u)*u_dot;           % Calculate Velocity
        yp_dot=r_dot*sin(u)+r*cos(u)*u_dot;           % Calculate Velocity

        omg = omg0(index) + (odot(index) - odote)*t - odote*toe(index);

        omg_dot=odot(index) - odote;                  % Calculate Velocity

        XS = xp*cos(omg) - yp*cos(i)*sin(omg);
        YS = xp*sin(omg) + yp*cos(i)*cos(omg);
        ZS = yp*sin(i);

        VXS = xp_dot*cos(omg)-yp_dot*cos(i)*sin(omg)...  % Calculate Velocity
            +yp*sin(i)*sin(omg)*i_dot-YS*omg_dot;
        VYS = xp_dot*sin(omg)+yp_dot*cos(i)*cos(omg)...  % Calculate Velocity
            -yp*sin(i)*i_dot*cos(omg)+XS*omg_dot;
        VZS = yp_dot*sin(i)+yp*cos(i)*i_dot;             % Calculate Velocity

        % compute the range
        R = [XS,YS,ZS] - XYZ_station(:,ti)';
        R = sqrt(sum(R.^2));
        tau=R/lightspeed;
        % earth rotation correction 
        phi=-odote*tau;
        corr=[cos(phi),-sin(phi);sin(phi),cos(phi)]*[XS;YS];
        XS=corr(1);
        YS=corr(2);
        % light Travel time correction
        XS = XS - VXS*tau;
        YS = YS - VYS*tau;
        ZS = ZS - VZS*tau;
        Xs(ii,ti) = XS;
        Ys(ii,ti) = YS;
        Zs(ii,ti) = ZS;
        Xv(ii,ti) = VXS;
        Yv(ii,ti) = VYS;
        Zv(ii,ti) = VZS;
    end
end

function E = kepOrb2E(M,e)
% Inputs:  - mean anomaly in radians
%          - eccentricity
% Output: Eccentric anomaly

if (-pi < M < 0) || (M > pi)
    E = M - e;
else
    E = M + e;
end

check = 1;

while check > 10e-10
    E_new = (E + (M - E + e * sin(E))/(1 - e * cos(E)));
    check = abs(E_new - E);
    E = E_new;
end
end

end