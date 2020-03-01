clear all
clc

addpath geom3d-2019.09.26/geom3d

% Part 1: define and generate the intersection

% Define geometry;
g.block.size.x = 250; % b
g.block.size.y = 200; % b
g.block.noBlocks.x = 3; % m_x
g.block.noBlocks.y = 2; % m_y, this is along the driving direction


% HEIGHT
g.build.meanHeight = 50; %log(20); % log(25); % for simulation: log(10,20,25)
g.build.stdDevHeight = 10; %log(3); %log(1.35);


% WIDTH
g.build.meanWidth = 50;
g.build.stdDevWidth = 0; %5

g.build.depth = 10;
g.build.qMax = 0; % max depth variation from street, uniform from 0 to this

g.road.width = 20; % \Delta b

disp('Generating urban canyon...')
last = 0;
blockCorners.x = NaN*ones(3,1);
for b = 1:1:g.block.noBlocks.x
    blockCorners.x(b) = last + g.road.width + g.block.size.x;
    last = blockCorners.x(b);
end
blockCorners.x = blockCorners.x - g.block.size.x;        % right corner to left corner

last = 0;
blockCorners.y = NaN*ones(3,1);
for b = 1:1:g.block.noBlocks.y
    blockCorners.y(b) = last + g.road.width + g.block.size.y;
    last = blockCorners.y(b);
end
blockCorners.y = blockCorners.y - g.block.size.y;            % blockCorners left lower point 

g.block.sizes = [g.block.size.x, g.block.size.y, 0]';     % ground

% loop for every block
P = zeros(0,1);
for bx = 1:1:g.block.noBlocks.x
    xFirst = bx == 1;       % if the block is the last one in x/y
    xLast = bx == g.block.noBlocks.x;
    xMiddle = ~xFirst && ~xLast;
    
    for by = 1:1:g.block.noBlocks.y
        
        if bx == 0 && by == 0
            % Skip middle building!!
        else
        
        yFirst = by == 1;
        yLast = by == g.block.noBlocks.y;
        yMiddle = ~yFirst && ~yLast;
        
        block.startPoint = [-95, -80, 0]' + [blockCorners.x(bx), blockCorners.y(by), 0]';
        
        % Decide which sides to skip, to avoid excessive amount of walls
        sideSkips(1).toSkip = [];
        sideSkips(2).toSkip = [];
        
        if yFirst
            sideSkips(1).toSkip = 1;
        elseif  yLast
            sideSkips(1).toSkip = 2;
        end
        
        if xFirst
            sideSkips(2).toSkip = 1;
        elseif xLast
            sideSkips(2).toSkip = 2;
        end
        
        
        for currDim = 1:2 % increasing x or y? -> decide which dimension to do
            otherDim = setdiff(1:2, currDim);
            currSkip =  sideSkips(currDim).toSkip;
            for currSide = setdiff(1:2,currSkip) % which side? -> decide which sides to do
                otherSide = setdiff(1:2, currSide);
                
                if currSide == 1
                    build.startPoint = [0, 0, 0]';
                elseif currSide == 2 && currDim == 1
                    build.startPoint = [0, g.block.sizes(2)-1*g.build.depth, 0]'; % y
                elseif currSide == 2 && currDim == 2
                    build.startPoint = [g.block.sizes(1)-1*g.build.depth, 0, 0]'; % x
                end
                
                build.startPoint = build.startPoint + block.startPoint;
                endCorner = 0;
                startCorner = 1;
                % Loop to fill with buildings
                while build.startPoint(currDim) < g.block.sizes(currDim) +  block.startPoint(currDim)
                    
                    
                    if g.build.stdDevWidth == 0
                        newWidth = g.build.meanWidth;
                    else
                        newWidth =  random('Rician', g.build.meanWidth, g.build.stdDevWidth);
                    end
                    
                    if newWidth > g.block.sizes(currDim) - build.startPoint(currDim) +  block.startPoint(currDim)
                        endCorner = 1;
                        newWidth = g.block.sizes(currDim) - build.startPoint(currDim) +  block.startPoint(currDim);
                    end
                    build.size = NaN*ones(3,1);
                    build.size(currDim) = newWidth;
                    
                    if g.build.stdDevHeight == 0
                        newHeight = g.build.meanHeight;
                    else
                        newHeight = random('Rician', g.build.meanHeight, g.build.stdDevHeight);
                    end
                    build.size(3) = newHeight;
                    
                    build.size(otherDim) = g.build.depth;
                    
                    build.upper = build.startPoint + build.size;
                    
                    
                    depthOffset_pos = 1*g.build.qMax*rand; % Uniform distribution
                    buildStartPointTemp = build.startPoint;
                    if currDim == 1 && currSide == 1
                        buildStartPointTemp(2) = buildStartPointTemp(2) + depthOffset_pos;
                    elseif currDim == 1 && currSide == 2
                        build.upper(2) = build.upper(2) - depthOffset_pos;
                    elseif currDim == 2 && currSide == 1
                        buildStartPointTemp(1) = buildStartPointTemp(1) + depthOffset_pos;
                    elseif currDim == 2 && currSide == 2
                        build.upper(1) = build.upper(1) - depthOffset_pos;
                    end
                    %                     if currSide == 1
                    %                         build.upper(otherDim) = build.upper(otherDim) + depthOffset;
                    %                     elseif currSide == 2
                    %                         build.startPoint(otherDim) = build.startPoint(otherDim) + depthOffset;
                    %                     end
                    
                    newPlanes = UCMP_generateBlock(buildStartPointTemp, build.upper, 'build');
                    
                    disp(build.upper);
                    disp(buildStartPointTemp);
                                        
                    plot3(buildStartPointTemp(1),buildStartPointTemp(2), buildStartPointTemp(3),'*r')
                    hold on 
                    plot3(build.upper(1),build.upper(2),build.upper(3), 'og')
                    hold on
                    grid on
                    
                    P = [P, newPlanes];
                    
                    newLims = [0, 0, 0]';
                    newLims(currDim) = build.size(currDim);
                    build.startPoint = build.startPoint + newLims;
                    startCorner = 0;
                end
            end
        end
        
        end
        
    end
end

upperPlotLim = max([blockCorners.x(end) + g.block.size.x + g.road.width...
    blockCorners.y(end) + g.block.size.y + g.road.width]);


noGeoPlanes = length(P)
disp('Urban canyon generated!')

fontSize = 40
% Attribute one color to each received PRN
possibleCols = {'r', 'g', 'b', 'm', 'c', 'k', 'y', [255,192,203]/255, [0.5 0.5 0.5], [255, 165, 0]/255, [0.1 0.1 0.1]};

blockColors.car = 0.8*ones(1,3);
blockColors.ground = 0.5*ones(1,3);
blockColors.build = 0.3*ones(1,3);

close all
figure('units', 'normalized', 'outerposition', [0 0 1 1])
box on
hold on
% Only last position of car will be plotted
for w = 1:1:length(P)
    poly_rectangle(P(w).points(:,1), P(w).points(:,2),...
        P(w).points(:,3), P(w).points(:,4),  blockColors.(P(w).objectType))     % fill3-----plot
    
    % Add plane ID number, for debuggung purposes (optional)
    %     text(P(w).lims(1,1) + (P(w).lims(1,2) -  P(w).lims(1,1))/2, ...
    %         P(w).lims(2,1) + (P(w).lims(2,2) -  P(w).lims(2,1))/2, ...
    %         P(w).lims(3,1) + (P(w).lims(3,2) -  P(w).lims(3,1))/2, ...
    %         num2str(w), 'FontSize', 15)
end

view(-10, 43)
xlabel('East (m)' , 'FontSize', fontSize, 'Interpreter', 'latex')
ylabel('North (m)', 'FontSize', fontSize, 'Interpreter', 'latex')
zlabel('Up (m)', 'FontSize', fontSize, 'Interpreter', 'latex')

zoomPoint = [0, 0, 0]'

% Option to zoom on car initial position:
zoomSize = 300; % m
%centerpoint = a.pos_init
xlim([zoomPoint(1)-zoomSize/3, zoomPoint(1)+zoomSize*4])
ylim([zoomPoint(2)-zoomSize, zoomPoint(2)+zoomSize])
zlim([0 zoomSize/2])

disp('Town square generated!')
%%

% Part 2: Import vehicle trajectory (sampled every 50 ms, only need every second)
spheroid = referenceEllipsoid('WGS84');
UMTfolder = './umt_obs'
% list of vehicles to use in simulation
%vehiclesToImport = [0:1:1] % to use all 15 vehicles! 
vehiclesToImport = [0];

vehColors = {'g', 'b', 'r', 'm', 'y'}
clear vehTraj
noV = length(vehiclesToImport)
for v = 1:1:noV
    
    % Read UMT files
    fileName = ['v' num2str(vehiclesToImport(v))];
    path = [UMTfolder '/' fileName '.umt'];
    
    parsedVehicle = parseUMT(path);
    vehTraj(v).t = parsedVehicle(1,:);
    vehTraj(v).ECEF = parsedVehicle(2:4,:);
    lla_linPoint = ecef2lla(vehTraj(1).ECEF(:,1)')'; % first vehicle decides where to lineraize
    
    vehTraj(v).LLA = ecef2lla(vehTraj(v).ECEF')';
    [xEast, yNorth, zUp] = ecef2enu(...
        vehTraj(v).ECEF(1,:), vehTraj(v).ECEF(2,:), vehTraj(v).ECEF(3,:),...
        lla_linPoint(1), lla_linPoint(2), lla_linPoint(3), spheroid);
    
    vehTraj(v).ENU = [xEast; yNorth; zUp];

end
disp('Vehicles UMT position data imported!')

%%
% (Optional: pre- plot the vehicle trajectories, to see that don't drive into
% the buildings
for v = 1:1:noV
    hold on
    plot3(vehTraj(v).ENU(1,:), vehTraj(v).ENU(2,:), ...
        vehTraj(v).ENU(3,:), ['-o' vehColors{v}])
end

%%

% [obs] = readRinexObsSprient("./umt_obs/rinex-obs_V1_A1-land vehicle.txt");
% car = CarModel("./umt_obs/stra_heading0.umt", obs);    % umt

% car = readtable('./umt_obs/car0.csv');
% car_maxTRow = height(car);

CLT = readtable('sat_data_V1A1.csv');                  % satellite log file 
%CLT = readtable('./umt_obs/sat_car0.csv');
sat_maxTRow = height(CLT);                            % Number of Lines in Log-Table

 
stp = 17;     % uniform linear motion 
lat0 = lla_linPoint(1);
lon0 = lla_linPoint(2);
h0 = lla_linPoint(3);
refl = 6;

%%
for i= 1:sat_maxTRow
                                  %%%%
    mint = floor((CLT.Time_ms(i))/60000);
    sec = mod((CLT.Time_ms(i))/1000, 60);
    [sow(i), tow(i)] = Date2GPSTime(2019,07,03,12, mint, sec);  %%%%%%%%%%%%%%%%%%% hour action time 
end

sate_log(:,1) = sow;
sate_log(:,2) = tow;
sate_log(:,3) = CLT.Sat_PRN;
sate_log(:,4) = CLT.Sat_Pos_X;
sate_log(:,5) = CLT.Sat_Pos_Y;
sate_log(:,6) = CLT.Sat_Pos_Z;
sate_log(:,7) = CLT.Azimuth;
sate_log(:,8) = CLT.Elevation;

% car_log(:,1) = car.TOW_ms/1000;
% car_log(:,2) = car.Pos_X;
% car_log(:,3) = car.Pos_Y;
% car_log(:,4) = car.Pos_Z;


[sat_e,sat_n,sat_u] = ecef2enu(sate_log(:,4),sate_log(:,5),sate_log(:,6),lat0,lon0,h0,wgs84Ellipsoid);   % satellite coordinate in ENU
sate_log(:,4:6) = [sat_e,sat_n,sat_u]; 

startp = find(sate_log(:,2) == car(stp,2));
sate_log = sate_log(startp(1):end,:);   % truncate
len_sate_log = length(sate_log);

%% plane building

H_wall = 3;  % Height of Wall in m
D_wall_r = 8;  % Distance to (East) right Wall in m
D_wall_l = 12; % Distance to (West) left Wand in m
L_cross = D_wall_r + D_wall_l;
len_wall = 500;

[xEast,yNorth,zUp] = ecef2enu(car(stp:end, 3),car(stp:end, 4),car(stp:end, 5),lat0,lon0,h0,wgs84Ellipsoid);

uni_motion = [car(stp:end,1:2), xEast, yNorth, zUp];    %%%%% car_matrix
% uni_motion = round(uni_motion,3);
cnt = length(uni_motion)

if rem(cnt,2) == 0            % distance of cross is 20m
    cross_st = cnt/2; 
else
    cross_st = (cnt-1)/2; 
end

cross1 = uni_motion(cross_st, 3:5);
cross2 = [cross1(1)+20,cross1(2),cross1(3)];

in_cross = find( uni_motion(:,3)>= cross1(1) & uni_motion(:,3) < cross2(1));       % only consider x direction 

%% Geometry
C2 = [cross1(1),cross1(2)+D_wall_l, cross1(3)];
C1 = [cross2(1),cross2(2)+D_wall_l, cross2(3)];
C4 = [cross2(1),cross2(2)-D_wall_r, cross2(3)];
C3 = [cross1(1),cross1(2)-D_wall_r, cross1(3)];

D1 = [C1(1)+len_wall,C1(2),C1(3)];
D4 = [C4(1)+len_wall,C4(2),C4(3)];
D2 = [C2(1)-len_wall,C2(2),C2(3)];
D3 = [C3(1)-len_wall,C3(2),C3(3)];

F1 = [C1(1),C1(2)+len_wall,C1(3)];
F2 = [C2(1),C2(2)+len_wall,C2(3)];
F3 = [C3(1),C3(2)-len_wall,C3(3)];
F4 = [C4(1),C4(2)-len_wall,C4(3)];

HC1 = [C1(1),C1(2),C1(3)+H_wall];
HC2 = [C2(1),C2(2),C2(3)+H_wall];
HC3 = [C3(1),C3(2),C3(3)+H_wall];
HC4 = [C4(1),C4(2),C4(3)+H_wall];

HD1 = [D1(1),D1(2),D1(3)+H_wall];
HD2 = [D2(1),D2(2),D2(3)+H_wall];
HD3 = [D3(1),D3(2),D3(3)+H_wall];
HD4 = [D4(1),D4(2),D4(3)+H_wall];

HF1 = [F1(1),F1(2),F1(3)+H_wall];
HF2 = [F2(1),F2(2),F2(3)+H_wall];
HF3 = [F3(1),F3(2),F3(3)+H_wall];
HF4 = [F4(1),F4(2),F4(3)+H_wall];

poly1_E = [C1; D1; HD1; HC1];    %11
poly1_N = [C1; F1; HF1; HC1];    %10
poly2_E = [C2; D2; HD2; HC2];    %21
poly2_N = [C2; F2; HF2; HC2];    %20
poly3_E = [C3; D3; HD3; HC3];    %31
poly3_N = [C3; F3; HF3; HC3];    %30
poly4_E = [C4; D4; HD4; HC4];    %41
poly4_N = [C4; F4; HF4; HC4];    %40

polys_all = zeros(4,3,8);
polys_all(:,:,1) = poly1_E;
polys_all(:,:,2) = poly1_N;
polys_all(:,:,3) = poly2_E;
polys_all(:,:,4) = poly2_N;
polys_all(:,:,5) = poly3_E;
polys_all(:,:,6) = poly3_N;
polys_all(:,:,7) = poly4_E;
polys_all(:,:,8) = poly4_N;

plane_w = createPlane(C2,C3,HC2);
plane_e = createPlane(C1,C4,HC1);
plane_n = createPlane(C1,C2,HC1);
plane_s = createPlane(C3,C4,HC3);

figure;
plot([C1(1),C2(1),C3(1),C4(1)],[C1(2),C2(2),C3(2),C4(2)],'*')
hold on
plot([D1(1),D2(1),D3(1),D4(1)],[D1(2),D2(2),D3(2),D4(2)],'o')
hold on
plot([F1(1),F2(1),F3(1),F4(1)],[F1(2),F2(2),F3(2),F4(2)],'.')


%% %%%%%LOS Detection%%%%%

LOSBlockedPoint = [];
LOSVisible = [];

for j = 1:cnt               %%%%%%%%%%%% starting point 
    key = find(sate_log(:,2) == uni_motion(j,2));
    car_pos = uni_motion(j, 3:5);
    lines = [];
    vis_sat = [];
    sat_prn = sate_log(key, 3);
    time = uni_motion(j,2);
    for i = 1:length(key)
        sat_pos = sate_log(key(i), 4:6);            %%%%%%%%%%%%%%%%%%%%  is the index right????????????????????
        sat_azi = sate_log(key(i), 7);
        sat_ele = sate_log(key(i), 8);
        

        L1 = createLine3d(car_pos, sat_pos);
        lines = [lines; L1];
        
        
%         if (abs(sat_azi)>0 & abs(sat_azi)<pi/2) | sat_azi> pi/2 
%             phi = abs(pi/2 - sat_azi);
%         else
%             phi = 3/2 * pi + sat_azi;

%         end
%         
%         L1 = createLine3d(car_pos,pi/2-sat_ele,phi); 
    end
            
    [PO1_E, inside1_E] = intersectLinePolygon3d(lines, poly1_E);
    PO1_E_PRN = [kron([time,11], ones(length(sat_prn(inside1_E)),1)), sat_prn(inside1_E), PO1_E(inside1_E,:)];    
    
    [PO1_N, inside1_N] = intersectLinePolygon3d(lines, poly1_N);
    PO1_N_PRN = [kron([time,10], ones(length(sat_prn(inside1_N)),1)), sat_prn(inside1_N), PO1_N(inside1_N,:)];
    
    [PO2_E, inside2_E] = intersectLinePolygon3d(lines, poly2_E);
    PO2_E_PRN = [kron([time,21], ones(length(sat_prn(inside2_E)),1)), sat_prn(inside2_E), PO2_E(inside2_E,:)];
    
    [PO2_N, inside2_N] = intersectLinePolygon3d(lines, poly2_N);
    PO2_N_PRN = [kron([time,20], ones(length(sat_prn(inside2_N)),1)), sat_prn(inside2_N), PO2_N(inside2_N,:)];
    
    [PO3_E, inside3_E] = intersectLinePolygon3d(lines, poly3_E);
    PO3_E_PRN = [kron([time,31], ones(length(sat_prn(inside3_E)),1)), sat_prn(inside3_E), PO3_E(inside3_E,:)];
    
    [PO3_N, inside3_N] = intersectLinePolygon3d(lines, poly3_N);
    PO3_N_PRN = [kron([time,30], ones(length(sat_prn(inside3_N)),1)), sat_prn(inside3_N), PO3_N(inside3_N,:)];
    
    [PO4_E, inside4_E] = intersectLinePolygon3d(lines, poly4_E);
    PO4_E_PRN = [kron([time,41], ones(length(sat_prn(inside4_E)),1)), sat_prn(inside4_E), PO4_E(inside4_E,:)];
    
    [PO4_N, inside4_N] = intersectLinePolygon3d(lines, poly4_N);
    PO4_N_PRN = [kron([time,40], ones(length(sat_prn(inside4_N)),1)), sat_prn(inside4_N), PO4_N(inside4_N,:)];

    
    inter = [PO1_E_PRN; PO1_N_PRN; PO2_E_PRN; PO2_N_PRN; PO3_E_PRN; PO3_N_PRN; PO4_E_PRN; PO4_N_PRN];
    
    if isempty(inter)
        vis_sat = find(sat_prn);
    else
        vis_sat = find(~ismember(sat_prn, inter(:,3)));
    end
    
    los = [kron(time, ones(length(vis_sat),1)), sat_prn(vis_sat), sate_log(key(vis_sat),7), sate_log(key(vis_sat),8)];
    
    LOSBlockedPoint = [LOSBlockedPoint; inter];    % timestamp, plane_code, prn, intersection point
    LOSVisible = [LOSVisible; los];
end

%% Reflection Detection

ReflectionPoint = [];
for j = 1:cnt
    
    key = find(sate_log(:,2) == uni_motion(j,2));
    car_pos = uni_motion(j, 3:5);
    lines = [];
    sat_prn = sate_log(key, 3);
    time = uni_motion(j,2);
    
    midpoint_n = projPointOnPlane(car_pos, plane_n);
    symPO_n = 2 * midpoint_n - car_pos;

    midpoint_e = projPointOnPlane(car_pos, plane_e);
    symPO_e = 2 * midpoint_e - car_pos;  

    midpoint_s = projPointOnPlane(car_pos, plane_s);
    symPO_s = 2 * midpoint_s - car_pos;
    
    midpoint_w = projPointOnPlane(car_pos, plane_w);
    symPO_w = 2 * midpoint_w - car_pos;
    
    if j < cross_st                         % left part of intersection 
        for i = 1:length(key)
            sat_pos = sate_log(key(i), 4:6);
            sat_azi = sate_log(key(i), 7);
            sat_ele = sate_log(key(i), 8);
            SVID = sate_log(key(i),3);
            
            if sat_azi > pi/2 || sat_azi < -pi/2    % car-> W to E, RHS, down part
                L1_point = symPO_n;               % mirror point1
                L1_poly1 = poly1_E;
                L1_poly2 = poly2_E;
                
                L2_point = symPO_e;               % mirror point2
                L2_interpoly = poly4_N;
                L2_noninter = poly3_N;
                L2_decide = poly4_E;
                otherpoly = poly3_E;
                                
            else
                L1_point = symPO_s;
                L1_poly1 = poly3_E;
                L1_poly2 = poly4_E;
                
                L2_point = symPO_e;
                L2_interpoly = poly1_N;
                L2_noninter = poly2_N;
                L2_decide = poly1_E; 
                otherpoly = poly2_E;
            end
            
            [refl, ReflectionPO, polycode, flag] = reflectionline(L1_point,L1_poly1, L1_poly2, L2_point,...
                                                        L2_interpoly, L2_noninter, L2_decide, otherpoly, sat_pos, cross2, L_cross, polys_all);
            
%             mirr_dist = norm(ReflectionPO - car_pos);
%             pathdelay = mirr_dist * sin(abs(sat_azi)) * cos(sat_ele);
            pathdelay = extrapath(car_pos, sat_pos, ReflectionPO);
                                                    
                                     
            if refl == 1
                reflect = [time, polycode, SVID, ReflectionPO, pathdelay];
                ReflectionPoint = [ReflectionPoint; reflect];                    %%%%  ReflectionPO means mirrorpoint
            end
            
                      
        end      
            
    elseif j> in_cross(end)               % right part of intersection 
        for i = 1:length(key)
            sat_pos = sate_log(key(i), 4:6);
            sat_azi = sate_log(key(i), 7);
            sat_ele = sate_log(key(i), 8);
            SVID = sate_log(key(i),3);
            
            if sat_azi > pi/2 || sat_azi < -pi/2
                L1_point = symPO_n;
                L1_poly1 = poly1_E;
                L1_poly2 = poly2_E;
                
                L2_point = symPO_w;
                L2_interpoly = poly3_N;
                L2_noninter = poly4_N;
                L2_decide = poly3_E;
                otherpoly = poly4_E;
                
            else
                L1_point = symPO_s;
                L1_poly1 = poly4_E;
                L1_poly2 = poly3_E;
                
                L2_point = symPO_w;
                L2_interpoly = poly2_N;
                L2_noninter = poly1_N;
                L2_decide = poly2_E;
                otherpoly = poly1_E;
            end
            
            [refl, ReflectionPO, polycode, flag] = reflectionline(L1_point,L1_poly1, L1_poly2, L2_point,...
                                                        L2_interpoly, L2_noninter, L2_decide, otherpoly, sat_pos, cross1, L_cross, polys_all);
        
%             mirr_dist = norm(ReflectionPO - car_pos);
%             pathdelay = mirr_dist * sin(abs(sat_azi)) * cos(sat_ele);                                                    
            pathdelay = extrapath(car_pos, sat_pos, ReflectionPO);                                      
                                                    
            if refl == 1
                reflect = [time, polycode, SVID, ReflectionPO, pathdelay];
                ReflectionPoint = [ReflectionPoint; reflect];
            end        

        end 
        
    else
        for i = 1:length(key)
            sat_pos = sate_log(key(i), 4:6);
            sat_azi = sate_log(key(i), 7);
            sat_ele = sate_log(key(i), 8);
            SVID = sate_log(key(i),3);
            
            if sat_azi > pi/2 || sat_azi < -pi/2
                L_N = createLine3d(symPO_n, sat_pos);
                PO_N1 = intersectLinePolygon3d(L_N, poly1_E);
                PO_N2 = intersectLinePolygon3d(L_N, poly2_E);
            
                L_W = createLine3d(symPO_w, sat_pos);
                PO_W3 = intersectLinePolygon3d(L_W, poly3_N);
                checkpoint_W3 = intersectLinePolygon3d(L_W, poly4_N); 
              
                L_E = createLine3d(symPO_e, sat_pos);
                PO_N4 = intersectLinePolygon3d(L_E, poly4_N);
                checkpoint_N4 = intersectLinePolygon3d(L_E, poly3_N);
                
                if ismember(0, isnan([PO_W3])) && all(isnan([checkpoint_W3])== 1) % intersectionpoint exists    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    refl = 1;
                    polycode = polygoncode(poly3_N, polys_all);
                    ReflectionPO = symPO_w;
                end
                
                if ismember(0, isnan([PO_N4])) && all(isnan([checkpoint_N4])== 1)  
                    refl = 1;
                    polycode = polygoncode(poly4_N, polys_all);
                    ReflectionPO = symPO_e;
                end                                                            %%%%%%%%%%%
                                
                
                if 1 && all(isnan([PO_N1,PO_N2])== 1)  % no intersection point
                    refl = 0;
                else
                    PO_check1 = intersectLinePolygon3d(L_N, poly4_E);
                    PO_check2 = intersectLinePolygon3d(L_N, poly4_N);
                    PO_check3 = intersectLinePolygon3d(L_N, poly3_N);
                    PO_check4 = intersectLinePolygon3d(L_N, poly3_E);
                    
                    if 1 && all(isnan([PO_check1, PO_check2, PO_check3, PO_check4])==1)
                        refl = 1;
                        if 1 && all(isnan(PO_N1)==1)    % no intersection point with poly1_E
                            polycode = polygoncode(poly2_E, polys_all);
                            ReflectionPO = symPO_n;
                        else
                            polycode = polygoncode(poly1_E, polys_all);
                            ReflectionPO = symPO_n;
                        end
                    else
                        refl = 0;
                    end
     
                end
                                    
               
            else
                L_S = createLine3d(symPO_s, sat_pos);
                PO_S4 = intersectLinePolygon3d(L_S, poly4_E);
                PO_S3 = intersectLinePolygon3d(L_S, poly3_E);
            
                L_W = createLine3d(symPO_w, sat_pos);
                PO_W2 = intersectLinePolygon3d(L_W, poly2_N);   
                checkpoint_W2 = intersectLinePolygon3d(L_W, poly1_N);   
              
                L_E = createLine3d(symPO_e, sat_pos);
                PO_N1 = intersectLinePolygon3d(L_E, poly1_N);
                checkpoint_N1 = intersectLinePolygon3d(L_E, poly2_N);
                
                if ismember(0, isnan([PO_W2])) && all(isnan([checkpoint_W2])== 1)
                    refl = 1;
                    polycode = polygoncode(poly2_N, polys_all);
                    ReflectionPO = symPO_w;
                end
                
                if ismember(0, isnan([PO_N1])) && all(isnan([checkpoint_N1])== 1)
                    refl = 1;
                    polycode = polygoncode(poly1_N, polys_all);
                    ReflectionPO = symPO_e;
                end
                                
                
                if 1 && all(isnan([PO_S4,PO_S3])== 1)  % no intersection point
                    refl = 0;
                else
                    PO_check1 = intersectLinePolygon3d(L_S, poly1_E);
                    PO_check2 = intersectLinePolygon3d(L_S, poly1_N);
                    PO_check3 = intersectLinePolygon3d(L_S, poly2_N);
                    PO_check4 = intersectLinePolygon3d(L_S, poly2_E);
                    
                    if 1 && all(isnan([PO_check1, PO_check2, PO_check3, PO_check4])==1) % NO INTERSECTION POINT
                        refl = 1;
                        if 1 && all(isnan(PO_S4)==1)    % NO INTERSECTION POINT WITH poly4_E
                            polycode = polygoncode(poly3_E, polys_all);
                            ReflectionPO = symPO_s;
                        else                            % 
                            polycode = polygoncode(poly4_E, polys_all);
                            ReflectionPO = symPO_s;
                        end
                    else
                        refl = 0;
                    end
     
                end
            end
            
%             mirr_dist = norm(ReflectionPO - car_pos);
%             pathdelay = mirr_dist * sin(abs(sat_azi)) * cos(sat_ele);        %%%%%%%%%%%%%%%%%%%%%% equation ??
            pathdelay = extrapath(car_pos, sat_pos, ReflectionPO);
            
            
            if refl == 1
                reflect = [time, polycode, SVID, ReflectionPO, pathdelay];
                ReflectionPoint = [ReflectionPoint; reflect];
            end          
                               
        end
    end
    
    
       
end



%%
% satellite visibility status 0:Blocked; 1:LOS reception; 2:LOS + NLOS
% reception, namely MP; 3:Only NLOS reception


SateStatus = [];
for j = 1:cnt
    
    time = uni_motion(j, 2);
    key1_los = LOSVisible(find(LOSVisible(:,1) == time), 2);             % LOS_prn 
    key2_nlos = ReflectionPoint(find(ReflectionPoint(:,1)== time), 3);   % NLOS_prn
    
    indx = find(sate_log(:,2) == uni_motion(j,2));   
    sate_id = sate_log(indx, 3);          % all the PRN 
    
    for m = 1:length(indx)
        sate_position = sate_log(indx(m), 4:6);
        sate_azimuth = sate_log(indx(m), 7);
        sate_elevation = sate_log(indx(m), 8);
        poly_id = 00;
        interpo = [0,0,0];
        delay = 0;
                
        if ismember(sate_id(m), key1_los)
            if ismember(sate_id(m), key2_nlos)
                status = 2;
                poly_id = ReflectionPoint(find(ReflectionPoint(:,1) == time & ReflectionPoint(:, 3)== sate_id(m)), 2);
                interpo = ReflectionPoint(find(ReflectionPoint(:,1) == time & ReflectionPoint(:, 3)== sate_id(m)), 4:6);
                delay = ReflectionPoint(find(ReflectionPoint(:,1) == time & ReflectionPoint(:, 3)== sate_id(m)), 7);
            else
                status = 1;
            end
        else
            if ismember(sate_id(m), key2_nlos)
                status = 3;
                poly_id = ReflectionPoint(find(ReflectionPoint(:,1) == time & ReflectionPoint(:, 3)== sate_id(m)), 2);
                interpo = ReflectionPoint(find(ReflectionPoint(:,1) == time & ReflectionPoint(:, 3)== sate_id(m)), 4:6);
                delay = ReflectionPoint(find(ReflectionPoint(:,1) == time & ReflectionPoint(:, 3)== sate_id(m)), 7);
            else
                status = 0;
                
            end
        end
        SateStatus = [SateStatus; [time, sate_id(m), status, sate_azimuth, sate_elevation, poly_id, interpo, sate_position, delay]];
    end
end


            
        
%% %% skyplot
temp = 2;    % 2  80 92
car_temp = uni_motion(temp, 3:5);
time_temp = uni_motion(temp, 2);
sat_temp = find(SateStatus(:,1) == time_temp);
labels = cellstr(num2str(SateStatus(sat_temp,2)));

BlockedByWall = [];
for a = 0:0.1:2*pi
    for e = 0:0.1:pi/2
        line_temp = createLine3d(car_temp, pi/2-e, -(a-pi/2));
        check1 = intersectLinePolygon3d(line_temp, poly4_E);
        check2 = intersectLinePolygon3d(line_temp, poly4_N);
        check3 = intersectLinePolygon3d(line_temp, poly3_N);
        check4 = intersectLinePolygon3d(line_temp, poly3_E);
        check5 = intersectLinePolygon3d(line_temp, poly2_E);
        check6 = intersectLinePolygon3d(line_temp, poly2_N);
        check7 = intersectLinePolygon3d(line_temp, poly1_N);
        check8 = intersectLinePolygon3d(line_temp, poly1_E);
        
        if ismember(0, isnan([check1 check2 check3 check4 check5 check6 check7 check8]))
            wall_b = [rad2deg(a), rad2deg(e)];
            BlockedByWall = [BlockedByWall; wall_b];
        end      
           
    end
end
    

[uni, uni_id] = unique(SateStatus(sat_temp,3));
leg = strings(length(uni),1);              % legend
for j = 1:length(uni)
    switch uni(j)                                     %0:Blocked; 1:LOS reception; 2:LOS + NLOS
        case 0                                        % reception, namely MP; 3:Only NLOS reception
            leg(j) = 'Blocked';
        case 1
            leg(j) = 'Only LOS reception';
        
        case 2
            leg(j) = 'LOS+NLOS reception';
        case 3
            leg(j) = 'Only NLOS reception';
    end
end
        
figure;
plots(20) = skyplot(BlockedByWall(:,1), BlockedByWall(:,2), '.k');
alpha(.5)
hold on

for i = 1:length(sat_temp)   
    switch SateStatus(sat_temp(i),3)
        case 0
            clr = '*r';
        case 1
            clr = '*b';
        case 2 
            clr = '*g';
        case 3
            clr = '*m';
    end
    
    plots(i) = skyplot(rad2deg(wrapTo2Pi(SateStatus(sat_temp(i), 4))), rad2deg(SateStatus(sat_temp(i), 5)), clr, labels(i));
    hold on
end

legend([plots(uni_id)],leg);   % plot the trajectory of satellites observed by Wettzell


%%
figure;view(3);axis equal; grid on; 
%axis([550 610 -30 30 0 20]);  % 92
axis([car_temp(1)-60 car_temp(1)+60 -30 30 0 20]);  % 80
%axis([car_temp(1)-30 car_temp(1)+60 -30 30 0 20]);  % 20
%axis([D2(1)-10 D1(1)+10 F4(2)-10 F1(2)+10 0 10]);
%axis([-50 50 -50 50 0 10]);

for i = 1:length(sat_temp)
    e = SateStatus(sat_temp(i),5);
    a = SateStatus(sat_temp(i),4);
    sv = SateStatus(sat_temp(i),10:12);
    line_temp = createLine3d(car_temp, pi/2-e, -(a-pi/2));
    inter_temp = SateStatus(sat_temp(i), 7:9);
    
    switch SateStatus(sat_temp(i),3)
        case 0
            continue
        case 1          
            drawLine3d(line_temp, 'b');
            hold on
        case 2
            drawLine3d(line_temp, 'b');
            hold on
            ref1_temp = createEdge3d(inter_temp, sv);
            drawEdge3d(ref1_temp);
            hold on
            ref2_temp = createEdge3d(inter_temp, car_temp);
            drawEdge3d(ref2_temp);
            hold on
        case 3
            ref1_temp = createEdge3d(inter_temp, sv);
            drawEdge3d(ref1_temp);
            hold on
            ref2_temp = createEdge3d(inter_temp, car_temp);
            drawEdge3d(ref2_temp);            
    end
    
end
            

fillPolygon3d(poly1_E, 'k','LineStyle','none');alpha(.5);
hold on
fillPolygon3d(poly1_N, 'k','LineStyle','none');alpha(.5);
hold on
fillPolygon3d(poly2_E, 'k','LineStyle','none');alpha(.5);
hold on
fillPolygon3d(poly2_N, 'k','LineStyle','none');alpha(.5);
hold on
fillPolygon3d(poly3_E, 'k','LineStyle','none');alpha(.5);
hold on
fillPolygon3d(poly3_N, 'k','LineStyle','none');alpha(.5);
hold on
fillPolygon3d(poly4_E, 'k','LineStyle','none');alpha(.5);
hold on
fillPolygon3d(poly4_N, 'k','LineStyle','none');alpha(.5);
hold off            

xlabel('E');
ylabel('N');
zlabel('U');
zlim([0 10]);
% xlim([-50 50]);
% ylim([-50 50]);


%%

fid = fopen('./umt_obs/mycoorcanyon.act','wt');

wt_t0 = SateStatus(1,1);

for i = 1:length(SateStatus)
    TIME_STAMP = actiontime(SateStatus(i,1), wt_t0);
    SAT_ID = SateStatus(i,2);
    SAT_OFF = num2str(SateStatus(i, 13));
    REF_LOSS = '10,';
    
    switch(SateStatus(i,3))
        case 0           % 0:Blocked;   
            fprintf(fid, [TIME_STAMP ',SWITCH_SAT,V1_A1_VT_GPS,'  num2str(SAT_ID) ',2,']);
            fprintf(fid, '\n');
        case 1              % 1:LOS reception;
            fprintf(fid, [TIME_STAMP ',SWITCH_SAT,V1_A1_VT_GPS,'  num2str(SAT_ID) ',0,']);
            fprintf(fid, '\n');
        case 2              % 2:LOS + NLOS reception
            fprintf(fid, [TIME_STAMP ',SWITCH_SAT,V1_A1_VT_GPS,' num2str(SAT_ID) ',3,1,0,1,2,<41>0,-1,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,,<41>1,-1,' REF_LOSS SAT_OFF ',0.5,0.1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,,,0,0']);
            fprintf(fid, '\n');     
        case 3              % 3:only NLOS reception
            fprintf(fid, [TIME_STAMP ',SWITCH_SAT,V1_A1_VT_GPS,' num2str(SAT_ID) ',3,1,1,1,2,<41>0,-1,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,,<41>1,-1,' REF_LOSS SAT_OFF ',0.5,0.1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,,,0,0']);
            fprintf(fid, '\n');
    end
end


fclose(fid);
SateStatus_duration = SateStatus(end,1) - SateStatus(1,1)


fid = fopen("./umt_obs/stra_heading0.umt",'r');
str = textscan(fid,'%s','Delimiter','\n');
fclose(fid);


car_t0 = car(1,2);   
car_stp = car(stp,2);


for i = 1:size(str{1})
    time_str = str2num(str{1}{i}(1:9)); 
    
    if mod(time_str,1) == 0
        umt_t0 = time_str;               %%% make the first value coincide 
        fixed_diff = car_t0 - umt_t0;
        umt_stp = (car_stp - fixed_diff) *1000;        
        break
    end
end

for i = 1:size(str{1})
    compr_str = str2num(str{1}{i}(1:9)) *1000;
    if compr_str == umt_stp
        trunc = i;
        break
    end
end


% Extract from line stp
str2 = str{1}(trunc:end);
umt_duration = str2num(str2{end}(1:9)) - str2num(str2{1}(1:9))


% Save as a text file
fid2 = fopen('./umt_obs/mycoorcanyon.umt','w');
fprintf(fid2,'%s\n', str2{:});
fclose(fid2);




% 
% 
% figure;view(3);axis equal; grid on; 
% axis([D2(1) D1(1) F4(2) F1(2) -100 100])
% 
% fillPolygon3d(poly1_E, 'k','LineStyle','none');alpha(.5);
% hold on
% fillPolygon3d(poly1_N, 'k','LineStyle','none');alpha(.5);
% hold on
% fillPolygon3d(poly2_E, 'k','LineStyle','none');alpha(.5);
% hold on
% fillPolygon3d(poly2_N, 'k','LineStyle','none');alpha(.5);
% hold on
% fillPolygon3d(poly3_E, 'k','LineStyle','none');alpha(.5);
% hold on
% fillPolygon3d(poly3_N, 'k','LineStyle','none');alpha(.5);
% hold on
% fillPolygon3d(poly4_E, 'k','LineStyle','none');alpha(.5);
% hold on
% fillPolygon3d(poly4_N, 'k','LineStyle','none');alpha(.5);
% hold on






% 
% c = 23      
% a = sate_log(c,7);
% a1 = wrapTo2Pi(a)
% e = sate_log(c,8)
% line1 = createLine3d([0,0,0],sate_log(c,4:6));
% line2 = createLine3d([0,0,0], pi/2-e, -(a1-pi/2))  %wrapTo2Pi(a-pi/2)
% 
% if (abs(a)>0 & abs(a)<pi/2) | a> pi/2 
%     phi = abs(pi/2 - a);
% else
%     phi = 3/2 * pi + a;
% end
% 
% line3 = createLine3d([0,0,0],pi/2-e,phi)
% 
% figure;view(3);axis equal; grid on; axis([-100 100 -100 100 -100 100])
% drawLine3d(line1,'k')
% hold on
% drawLine3d(line2, 'b') 
% % hold on
% % drawLine3d(line3, 'g')     
% xlabel('x')
% ylabel('y')
% zlabel('z')


%         if (abs(sat_azi)>0 & abs(sat_azi)<pi/2) | sat_azi> pi/2 
%             phi = abs(pi/2 - sat_azi);
%         else
%             phi = 3/2 * pi + sat_azi;
%         end
%         
%         L1 = createLine3d(car_pos,pi/2-sat_ele,phi); 

                

 %%                     
%         [sat_e,sat_n,sat_u] = ecef2enu(sat_pos(1),sat_pos(2),sat_pos(3),lat0,lon0,h0,wgs84Ellipsoid);
%         sat_enu = [sat_e,sat_n,sat_u]
%         L2 = createLine3d(car_pos,sat_enu);drawLine3d(L2,'b');
%                      
%         L1 = createLine3d([0,0,0],atan(sqrt(2)),pi*3/4)
%         L2 = createLine3d([0,0,0],[-1,1,1])



%         figure;view(3);axis equal; grid on; axis([D2(1) D1(1) F4(2) F1(2) -100 100])
%         drawLine3d(L1,'k')
%         hold on       
%         drawLine3d(L2,'b')
%         light;
%         xlabel('x');
%         ylabel('y');
%         zlabel('z');
%       
%         plane1 = createPlane(C1,D1,HC1);
%         plane2 = createPlane(C2,D2,HC2);
%         L1 = createLine3d(C1,D1);
%         L2 = createLine3d(C2,D2);
% 
%         axis([D2(1) D1(1) F4(2) F1(2) -100 100]);
%         drawPlane3d(plane1)
%         drawPlane3d(plane2,'b')
%         drawLine3d(L1,'b')
%         drawLine3d(L2,'k')
%         set(gcf, 'renderer', 'zbuffer');
%         xlabel('x');
%         ylabel('y');
%         zlabel('z');
%         
% 
%         
%         
%         PO = intersectLinePlane(L, plane1)
%         h = drawPlane3d(plane1)
%         
%         figure; drawLine3d(L2)
%         
%     
%     
%     

