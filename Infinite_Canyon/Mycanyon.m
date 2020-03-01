
clear;
close all;

[obs] = readRinexObsSprient("./umt_obs/rinex-obs_V1_A1-land vehicle.txt");
car = CarModel("./umt_obs/stra_heading0.umt", obs);    % umt

CLT = readtable('sat_data_V1A1.csv');                  % satellite log file 
maxTRow = height(CLT);                            % Number of Lines in Log-Table


PRN = unique(CLT.Sat_ID)
timestamp = unique(CLT.Time_ms)
stp = 15     % uniform linear motion 
lat0 = 48.726247
lon0 = 9.067573
h0 = 0
%%
figure;
plot3(car(:,3), car(:,4), car(:,5))
figure;
plot(car(stp:end,3), car(stp:end,4))

hold on 
figure;
plot3(sate_log(:,4), sate_log(:,5), sate_log(:,6),'.')
hold on
plot3(car(:,3), car(:,4), car(:,5))
grid on
%% plane building

[xEast,yNorth,zUp] = ecef2enu(x,y,z,lat0,lon0,h0,wgs84)


A = car(stp,3:5);
B = car(end,3:5);
dif = B - A;
e = dif/norm(dif); 
n1 = [-e(2), e(1)];
n2 = [e(2), -e(1)];
MD = (A + B)/2;

%%
for i= 1:maxTRow
    mint = floor((CLT.Time_ms(i))/60000);
    sec = mod((CLT.Time_ms(i))/1000, 60);
    [sow(i), tow(i)] = Date2GPSTime(2019,07,03,12, mint, sec);
end

sate_log(:,1) = tow';
sate_log(:,2) = sow';
sate_log(:,3) = CLT.Sat_PRN;
sate_log(:,4) = CLT.Sat_Pos_X;
sate_log(:,5) = CLT.Sat_Pos_Y;
sate_log(:,6) = CLT.Sat_Pos_Z;
sate_log(:,7) = CLT.Azimuth;
sate_log(:,8) = CLT.Elevation;

%%
H_wall = 10;  % Height of Wall in m
D_wall_r = 5;  % Distance to (East) right Wall in m
D_wall_l = 10;  % Distance to (West) left Wand in m

[xEast,yNorth,zUp] = ecef2enu(car(stp:end, 3),car(stp:end, 4),car(stp:end, 5),lat0,lon0,h0,wgs84Ellipsoid);

uni_motion = [car(stp:end,1:2), xEast, yNorth, zUp];
cnt = length(uni_motion);

if rem(cnt,2) == 0            % distance of cross is 15m
    cross_st = cnt/2 
else
    cross_st = (cnt-1)/2 
end

cross1 = [xEast(cross_st), yNorth(cross_st), zUp(cross_st)];
cross2 = [xEast(cross_st) + 20, yNorth(cross_st), zUp(cross_st)];

dis1 = round([xEast, yNorth, zUp] - kron(cross1, ones(cnt,1)),2);    %%%%%%%% round or not?
dis2 = round([xEast, yNorth, zUp] - kron(cross2, ones(cnt,1)),2);

crspoint = find(dis1(:,1)>=0 & dis1(:,1)<=20);
cross_ed = min(find(dis1(:,1)>20));
% figure;
% plot3(xEast, round(yNorth,2), round(zUp,2))
% figure;
% plot(xEast,round(yNorth,2))

azi = sate_log(:,7);
ele = sate_log(:,8);


for i = 1:cross_st
    key = find(sate_log(:,2) == uni_motion(i,2));
    for j = 1:length(key)                % j is for satellite related; i is car related
        if azi(key(j)) > 0   % east
            min_El_LOS = atan(H_wall/(D_wall_r/sin(azi(key(j)))));   % min. Elevation-Angle for LOS Reception 
            min_azi_LOS = atan(-dis1(i)/D_wall_r) + pi/2 
            max_azi_LOS = atan(-dis2(i)/D_wall_r) + pi/2
            min_ele_azi = atan(H_wall/abs(dis2(i)/cos(azi(key(j)))))
            
            if ele(key(j)) > min_EL_LOS || (azi(key(j)) > min_azi_LOS && azi(key(j)) < max_azi_LOS && ele(key(j)) > min_ele_azi)
                LOS = 1
            else
                LOS = 0
            end
            
            min_El_Refl = atan(sin(azi(key(j)))*H_wall/(2*D_wall_l+D_wall_r));   % min. Elevation Angle of 1st Order Reflection Reception
            max_El_Refl = atan(sin(azi(key(j)))*H_wall/D_wall_l);                % max. Elevation Angle of 1st Order Reflection Reception
            
            min_azi_Refl = atan(-dis1(i)/D_wall_l) + pi/2 
            max_azi_Refl = atan(-dis2(i)/D_wall_l) + pi/2
            
            azi1_refl = atan(D_wall_r/-dis2(i)) 
            azi2_refl = atan(D_wall_r/(-dis2(i)-20))
            
            if (ele(key(j)) > min_El_Refl) && (ele(key(j)) < max_El_Refl) && (azi(key(j))< min_azi_Refl) && (azi(key(j)) > max_azi_Refl) || ...
                    (ele(key(j)) < min_ele_azi && azi(key(i))< pi/2)
                
                Refl = 1
            else 
                Refl = 0
            end
            
        else   % west
            azi(key(i)) = abs(azi(key(i)))
            
            min_El_LOS = atan(H_wall/(D_wall_l/sin(azi(key(j)))));   % min. Elevation-Angle for LOS Reception 
            min_azi_LOS = atan(-dis1(i)/D_wall_l) + pi/2 
            max_azi_LOS = atan(-dis2(i)/D_wall_l) + pi/2
            min_ele_azi = atan(H_wall/abs(dis2(i)/cos(azi(key(j)))))
            
            if ele(key(j)) > min_EL_LOS || (azi(key(j)) > min_azi_LOS && azi(key(j)) < max_azi_LOS && ele(key(j)) > min_ele_azi)
                LOS = 1
            else
                LOS = 0
            end
            
            min_El_Refl = atan(sin(azi(key(j)))*H_wall/(2*D_wall_r+D_wall_l));   % min. Elevation Angle of 1st Order Reflection Reception
            max_El_Refl = atan(sin(azi(key(j)))*H_wall/D_wall_r);                % max. Elevation Angle of 1st Order Reflection Reception
            
            max_azi_Refl = -(atan(-dis1(i)/D_wall_r) + pi/2)       % -1.5
            min_azi_Refl = -(atan(-dis2(i)/D_wall_r) + pi/2)       %   -1.8
            
            
            if (ele(key(j)) > min_El_Refl) && (ele(key(j)) < max_El_Refl) && (azi(key(j))< min_azi_Refl) && (azi(key(j)) > max_azi_Refl) || ...
                    (ele(key(j)) < min_ele_azi && azi(key(i))< -pi/2)
                
                Refl = 1
            else 
                Refl = 0
            end
            
                
                

             




%%
for j = stp:size(car,1)               %%%%%%%%%%%% starting point 15s
    indx = find(sate_log(:,2) == car(j,2));
    car_pos = car(j, 3:5);
    for i = 1:length(indx)
        sat_pos = sate_log(i, 4:6);
        sat_azi = sate_log(i, 7);
        sat_ele = sate_log(i, 8);
        
        
    end
    
end




%%

lat1 = 0; 
lon1 = 0; 
alt1 = 1000*1000;
lat2 = 0; 
lon2 = 90; 
alt2 = 1000*1000;
[az, elev, r] = geodetic2aer(lat2,lon2,alt2,lat1,lon1,alt1,referenceEllipsoid('grs 80'));
Z = zeros(181,361);
R = georefpostings([-90 90],[-180 180], size(Z))
los2(Z,R,lat1,lon1,lat2,lon2,alt1,alt1); 
%%
Z = 500*peaks(100);
refvec = [1000 0 0];
[lat1, lon1, lat2, lon2] = deal(-0.027, 0.05, -0.093, 0.042);
los2(Z,refvec,lat1,lon1,lat2,lon2,100);

figure;
axesm('globe','geoid',earthRadius('meters'))
meshm(Z, refvec, size(Z), Z); axis tight
camposm(-10,-10,1e6); camupm(0,0)
demcmap('inc', Z, 1000); shading interp; camlight
[vis,visprofile,dist,h,lattrk,lontrk] = ... 
los2(Z,refvec,lat1,lon1,lat2,lon2,100);
plot3m(lattrk([1;end]),lontrk([1; end]),...
h([1; end])+[100; 0],'r','linewidth',2)
plotm(lattrk(~visprofile),lontrk(~visprofile),...
h(~visprofile),'r.','markersize',10)
plotm(lattrk(visprofile),lontrk(visprofile),...
h(visprofile),'g.','markersize',10)

%%
%[y z] = meshgrid(-1:0.1:1);
% syms x(y)
% x(y) = piecewise(y > 5, 5, y< -10, 5)
% fplot(x)

c = 0:0.1:7;    %height
b1 = 5:0.1:100;
b2 = -100:0.1:-10;
%b = [b1,b2];
[y1,z1] = meshgrid(b1, c);
[y2,z2] = meshgrid(b2, c);

x1_e = ones(size(y1))* 5;
x2_e = ones(size(y2))* 5;

x1_w = ones(size(y1))* -10;
x2_w = ones(size(y2))* -10;

% d1 = 5:0.1:100;
% d2 = -100:0.1:-10;
[xx1,zz1] = meshgrid(b1, c);
[xx2,zz2] = meshgrid(b2, c);

yy1_n = ones(size(xx1)) * 5;
yy2_n = ones(size(xx2)) * 5; 

yy1_s = ones(size(xx1)) * -10;
yy2_s = ones(size(xx2)) * -10;

figure();
%s = surf(x1_e, y1, z1,'FaceColor','b')
%s.EdgeColor = 'none';
surf(x1_e, y1, z1);
hold on
surf(x1_w, y1, z1);
hold on
surf(x2_e, y2, z2);
hold on
surf(x2_w, y2, z2);
hold on
surf(xx1, yy1_n, zz1)
hold on
surf(xx1, yy1_s, zz1)
hold on
surf(xx2, yy2_n, zz2)
hold on
surf(xx2, yy2_s, zz2)
xlim([-100,100])
ylim([-100,100])
zlim([0,20])


