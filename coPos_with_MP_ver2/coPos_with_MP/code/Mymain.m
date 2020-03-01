init

% Part 1: define and generate the town square!

% Define geometry;
g.block.size.x = 55; % b
g.block.size.y = 55; % b
g.block.noBlocks.x = 3; % m_x
g.block.noBlocks.y = 3; % m_y, this is along the driving direction


% HEIGHT
g.build.meanHeight = 50 %log(20); % log(25); % for simulation: log(10,20,25)
g.build.stdDevHeight = 10 %log(3); %log(1.35);


% WIDTH
g.build.meanWidth = 11;
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
blockCorners.x = blockCorners.x - g.block.size.x;

last = 0;
blockCorners.y = NaN*ones(3,1);
for b = 1:1:g.block.noBlocks.y
    blockCorners.y(b) = last + g.road.width + g.block.size.y;
    last = blockCorners.y(b);
end
blockCorners.y = blockCorners.y - g.block.size.y;

g.block.sizes = [g.block.size.x, g.block.size.y, 0]';

P = zeros(0,1);
for bx = 1:1:g.block.noBlocks.x
    xFirst = bx == 1;
    xLast = bx == g.block.noBlocks.x;
    xMiddle = ~xFirst && ~xLast;
    
    for by = 1:1:g.block.noBlocks.y
        
        if bx == 2 && by == 2
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
                    build.startPoint = [0, g.block.sizes(2)-1*g.build.depth, 0]';
                elseif currSide == 2 && currDim == 2
                    build.startPoint = [g.block.sizes(1)-1*g.build.depth, 0, 0]';
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

% Add the ground as reflection plane:
extraPlaneLims = 1e4;
groundStruct.normal = [0, 0, 1]';
groundStruct.lims = [-extraPlaneLims upperPlotLim+extraPlaneLims;...
    -extraPlaneLims upperPlotLim+extraPlaneLims;...
    -extraPlaneLims extraPlaneLims];
groundStruct.points = [-extraPlaneLims, -extraPlaneLims, 0 ;...
    upperPlotLim+extraPlaneLims, -extraPlaneLims, 0;...
    upperPlotLim+extraPlaneLims, upperPlotLim+extraPlaneLims, 0;...
    -extraPlaneLims, upperPlotLim+extraPlaneLims, 0]';
groundStruct.planeType = 'floor';
groundStruct.objectType = 'ground';
P(end+1) = groundStruct;

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
        P(w).points(:,3), P(w).points(:,4),  blockColors.(P(w).objectType))
    
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
xlim([zoomPoint(1)-zoomSize/2, zoomPoint(1)+zoomSize/2])
ylim([zoomPoint(2)-zoomSize/2, zoomPoint(2)+zoomSize/2])
zlim([0 zoomSize])

disp('Town square generated!')
%%

% Part 2: Import vehicle trajectory (sampled every 50 ms, only need every second)
spheroid = referenceEllipsoid('WGS84');
UMTfolder = '../../Infinite_Canyon/umt_obs/motion_V1.csv'
car = readtable(UMTfolder);
% list of vehicles to use in simulation
 vehiclesToImport = [1] % to use all 15 vehicles! 

vehColors = {'g', 'b', 'r'}
clear vehTraj
noV = length(vehiclesToImport)
for v = 1:1:noV
    
    % Read UMT files
%     fileName = ['v' num2str(vehiclesToImport(v))];
%     path = [UMTfolder '/' fileName '.umt'];
    
    %parsedVehicle = parseUMT(path);
    vehTraj(v).t = car.TOW_ms/1000;
    vehTraj(v).ECEF = [car.Pos_X, car.Pos_Y, car.Pos_Z]';
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
% Part 3: Compute reflections (linear algebra part)
% This takes around 2 minutes per vehicle, i.e., 30 minutes for all 15
% vehicles. The results are saved at the end, so one does not have to redo
% this step when after a Matlab restart

% Reflection settings, must be at start
sampTime = 1 % T, sampling rate, s % 1 Hz, like Rinex OBS files
sampTime_veh = 0.05
sampFreq = 1/sampTime;

% Car dimensions
carHeight = 1.5
carLength = 2%4.7
carWidth = 2%1.8

antAboveRoof = 0.01 % antenna heigth, delta 0.05

antHeight = carHeight +  antAboveRoof % vehTraj_us(3,1);



startTime = 86400/sampFreq
% Generate satellite trajectories given this start time

% Option 1 for satellite positions: generate new satellite positions
% [satPos_ECEF_X, satPos_ECEF_Y, satPos_ECEF_Z,...
%     satVel_ECEF_X, satVel_ECEF_Y, satVel_ECEF_Z] =...
%     UCMP_generateSatTraj_ECEF(startTime, 1/sampFreq, simDuration);


% Option 2 for satellite positions: read rinex and convert to mat
%filename_nav = [ '../data/square/20190703.nav'];
%ephemeris = readRinexNav(filename_nav);

% Option 3 for satellite positions: load .mat file (currently the best)
filename_nav = [ '../data/square/03072019.mat'];
ephemeris_load = load(filename_nav);
ephemeris = ephemeris_load.ephemeris;
    
clear res
close all
logicProblem = 0;
stepCompTime = [];
disp('Starting computation of reflections...')
for v = 1:1:noV % Loop over all vehicles
    
    % Undersample the vehicle trajectory
    [~, noSamp] = size(vehTraj(v).ENU)
    vehTraj_us = vehTraj(v).ENU(:,1:sampTime/sampTime_veh:end) + [0; 0; antHeight];
    vehTime_us = vehTraj(v).t(:,1:sampTime/sampTime_veh:end);

    antVel_init_ENU = vehTraj_us(:,2) - vehTraj_us(:,1);
    
    antLat_init_deg = vehTraj(v).LLA(1,1); % degrees, lambda
    antLon_init_deg = vehTraj(v).LLA(2,1); % degrees, phi
    
    % Computed settings:
    
    antPos_init_ENU = vehTraj_us(:,1) 
    antPos_init_LLA = vehTraj(v).LLA(:,1)
    
    % Simulate until the middle of the last block
    startTime_proc = 1 % We can stop at 50, since 1=51 in terms of constellation. percentage of day, 0 to 100, 50
    %simDistance = (g.block.noBlocks.y - 1)*g.block.size.y + (g.block.noBlocks.y - 1)*g.road.width;
    
    simDuration = vehTime_us(end) - vehTime_us(1) + 1  % seconds
    
    XYZ_station = vehTraj(v).ECEF;
 [satPos_ECEF_X, satPos_ECEF_Y, satPos_ECEF_Z,...
     satVel_ECEF_X, satVel_ECEF_Y, satVel_ECEF_Z] = ...
    UCMP_generateSatTraj_ECEF_Yan(startTime, sampTime, simDuration, ...
    ephemeris, XYZ_station);
    
    % Option to rotate the entire scenario
    % g.rot = 0 % of rot = 0, East = x, North = y, Up = z
    
    no_samps = length(vehTime_us)
    totalSteps = no_samps * noV
    
    % Satellite trajectory parameters
    s.allPRN = 1:1:31;
    
    antPos_prev_ENU = antPos_init_ENU;
    

    for samp = 1:1:no_samps
        tic
        
        currTime = startTime + samp*sampTime;
        res(v, samp).currTime = currTime;
        
        % Update antenna position
        %  antHead_rad =  1
        %
        %     velRotMatr = inv([cos(antHead_rad) -1*sin(antHead_rad)  0 ;
        %         sin(antHead_rad) cos(antHead_rad) 0;
        %         0 0 1]);
        
        %
        antPos_ENU = vehTraj_us(:,samp);
        %     antPos_prev_ENU = antPos_ENU;
        %     res(constRep, samp).antHead_rad = antHead_rad;
        %     res(constRep, samp).velRotMatr = velRotMatr;
        res(v, samp).antPos_ENU = antPos_ENU;
        res(v, samp).antPos_prev_ENU = antPos_prev_ENU;
        antVel_ENU = antPos_ENU -  antPos_prev_ENU;
        res(v, samp).antVel_ENU = antVel_ENU;
        
        % Antenna position in longitude, latitude_altitude
        [antPos_LLA(1), antPos_LLA(2), antPos_LLA(3)] = enu2geodetic(antPos_ENU(1), antPos_ENU(2), antPos_ENU(3),...
            antPos_init_LLA(1), antPos_init_LLA(2), antPos_init_LLA(3), spheroid);
        res(v, samp).antPos_LLA = antPos_LLA';
        
        % Must also rotate the car (or keep length of car = width of car)!
        carLower = antPos_ENU - 1*[carWidth/2 carLength/2 antHeight]'; % Car is slighly lower than antenna
        carUpper = antPos_ENU + [carWidth/2 carLength/2 carHeight-antHeight]';
        currCar = UCMP_generateBlock(carLower, carUpper, 'car'); % Last slot is always the car
        
        P(noGeoPlanes+1:end) = [];
        P = [P, currCar];
        
        % Loop through satellites, here parameters must be logged depending on PRN!
        r = 0; % Counts reflections for each PRN
        for currPRN = s.allPRN
            % Get satellite position and transform to ENU
            satPos_ENU = NaN(3,1);
            satPos_ECEF = [satPos_ECEF_X(currPRN, samp); satPos_ECEF_Y(currPRN, samp); satPos_ECEF_Z(currPRN, samp)];
            [satPos_ENU(1), satPos_ENU(2), satPos_ENU(3)] = ecef2enu(satPos_ECEF(1), satPos_ECEF(2), satPos_ECEF(3), ...
                antPos_LLA(1), antPos_LLA(2), antPos_LLA(3), spheroid, 'degrees');
            
            satPos_LLA = ecef2lla(satPos_ECEF')';
            [satAzim_deg, satElev_deg, satSlantRange] = geodetic2aer(satPos_LLA(1), satPos_LLA(2), satPos_LLA(3), ...
                antPos_LLA(1), antPos_LLA(2), antPos_LLA(3), spheroid, 'degrees');
            satPos_AER = [satAzim_deg, satElev_deg, satSlantRange]'; % In degrees, degrees, metres!
            
            % Get satellite velocity and transform to ENU
            satVel_ECEF = [satVel_ECEF_X(currPRN, samp); ...
                satVel_ECEF_Y(currPRN, samp); satVel_ECEF_Z(currPRN, samp)];
            rotMat_ECEF_to_ENU =...
                [-sind(satElev_deg)                 cosd(satElev_deg)                   0;...
                -sind(satAzim_deg)*cosd(satElev_deg)    -sind(satAzim_deg)*sind(satElev_deg)    cosd(satAzim_deg);...
                cosd(satAzim_deg)*cosd(satElev_deg)     cosd(satAzim_deg)*sind(satElev_deg)     sind(satAzim_deg)];
            satVel_ENU = (rotMat_ECEF_to_ENU * satVel_ECEF);
            
            res(v, samp).satPos_ECEF(1:3, currPRN) = satPos_ECEF;
            res(v, samp).satPos_ENU(1:3, currPRN) = satPos_ENU;
            res(v, samp).satPos_AER(1:3, currPRN) = satPos_AER;
            res(v, samp).satVel_ECEF(1:3, currPRN) = satVel_ECEF;
            res(v, samp).satVel_ENU(1:3, currPRN) = satVel_ENU;
            
            % Check if line of sight is blocked:
            satLOS = 1;
            satArrival = antPos_ENU;
            for w = 1:1:length(P)
                % fisk
                crossObstacle = UCMP_planeLineIntersect(P(w).normal,...
                    P(w).points(:,1), P(w).lims, antPos_ENU, satPos_ENU);
                if sum(isnan(crossObstacle)) == 0 % Block detected, LOS is blocked!
                    satLOS = 0;
                    satArrival = crossObstacle;
                    break
                end
            end
            res(v, samp).satLOS(currPRN) = satLOS;
            res(v, samp).satArrival(1:1:3, currPRN) = satArrival;
            
            % Compute and add reflections for this time step
            for w1 = 1:1:length(P)
                % apa
                [reflPoint, mirrorPoint, projPoint] = UCMP_computeReflection(P(w1).normal,...
                    P(w1).points(:,1), P(w1).lims, antPos_ENU, satPos_ENU);
                
                if sum(isnan(reflPoint)) == 0  % A reflection has been found!
                    
                    % Check all other walls for crossings, to see if the
                    % reflection is blocked:
                    reflBlocked = 0;
                    for w2 = setdiff(1:1:length(P),w1)
                        crossObstacle_toBuild = UCMP_planeLineIntersect(P(w2).normal,...
                            P(w2).points(:,1), P(w2).lims, satPos_ENU, reflPoint);
                        crossObstacle_fromBuild = UCMP_planeLineIntersect(P(w2).normal,...
                            P(w2).points(:,1), P(w2).lims, reflPoint, antPos_ENU);
                        
                        % Reflection block detected, do not consider reflection!
                        if  (sum(isnan(crossObstacle_toBuild)) == 0) || (sum(isnan(crossObstacle_fromBuild)) == 0 ) %
                            reflBlocked = 1;
                            break
                        end
                    end
                    
                    if reflBlocked == 0 % Reflection found and it is not blocked, add to list
                        
                        r = r + 1; % Increase reflection count for this constRep and samp
                        res(v, samp).refls(r).prn = currPRN;
                        res(v, samp).refls(r).reflPoint = reflPoint;
                        res(v, samp).refls(r).mirrPoint = mirrorPoint;
                        res(v, samp).refls(r).projPoint = projPoint;
                        
                        res(v, samp).refls(r).planeType = P(w1).planeType;
                        res(v, samp).refls(r).objectType = P(w1).objectType;
                        
                        % Compute and check reflection angle
                        vec1 = satPos_ENU - reflPoint; % Order here is important!
                        vec2 = antPos_ENU - reflPoint;
                        planeNorm = P(w1).normal;
                        
                        inAngle_rad = pi/2 - acos( sum(vec1.* planeNorm) / (norm(vec1)*norm(planeNorm)) );
                        outAngle_rad = pi/2 - acos( sum(vec2.* planeNorm) / (norm(vec2)*norm(planeNorm)) );
                        
                        % inAngle and outAngle relative reflector plane must be equal and between 0 and 90!
                        if  inAngle_rad < 0 || inAngle_rad > pi/2 || abs(inAngle_rad - outAngle_rad) > 1e-2
                            disp('Reflection logic error!')
                            logicProblem = logicProblem + 1;
                        end
                        res(v, samp).refls(r).reflAngle_rad = inAngle_rad;
                        
                    end
                end
            end
            
        end
        
        res(v, samp).noRefl = r;
        
        % Dynamically estimate simulation time
        stepCompTime = [stepCompTime toc];
        if mod(samp,10) == 0 % Output time estimation every 10th samle
            remainingSteps = no_samps*(noV - v) + no_samps - samp;
            remMins = mean(stepCompTime) * remainingSteps / 60
        end
        
    end
    
    disp('Linear algebra computed!')
    if logicProblem > 0
        disp('Check equations, something is wrong...')
    end
end
tempData_folder = '../data/square/temp_data/'
tempData_fileName = 'testsquare_reflections'
save([tempData_folder tempData_fileName '.mat'])
disp('Results saved!')
%%

% Reload saved results from reflections
init
tempData_folder = '../data/square/temp_data/'
tempData_fileName = 'testsquare_reflections'
clear D
D(1).data = load([tempData_folder tempData_fileName '.mat']);

%% Compute reflections (physics part)
rej_LHCP = 3; % k, dB 
fc = [const.GPS.L1.CA.f_carrier, const.GPS.L2.CL.f_carrier, const.GPS.L5.SOL.f_carrier];

no_f = length(fc);

for d = 1:1:1
    
    for v = 1:1:D(d).data.noV
        for samp = 1:1:D(d).data.no_samps
            currR = D(d).data.res(v, samp);
            
            for r = 1:1:currR.noRefl
                
                prn = currR.refls(r).prn;
                
                extraPath = norm(currR.antPos_ENU - currR.refls(r).reflPoint)...
                    + norm(currR.satPos_ENU(:,prn) - currR.refls(r).reflPoint)...
                    - norm(currR.antPos_ENU - currR.satPos_ENU(:,prn));
                D(d).data.res(v, samp).refls(r).Delta_d(1:no_f) = extraPath;
                
                % Compute 3D Doppler of direct signal ( should be in
                % the above part)
                f_LOS = UCMP_computeDoppler(currR.antPos_ENU, currR.antVel_ENU,...
                    currR.satPos_ENU(:, prn), currR.satVel_ENU(:, prn), fc, const);
                
                % Compute 3D Doppler of reflection
                f_toBuild = UCMP_computeDoppler(currR.refls(r).reflPoint, [0, 0, 0]',...
                    currR.satPos_ENU(:,prn), currR.satVel_ENU(:,prn), fc, const);
                
                f_refl = UCMP_computeDoppler(currR.antPos_ENU, currR.antVel_ENU,...
                    currR.refls(r).reflPoint, [0, 0, 0]', f_toBuild, const);
                
                dopplerDiff = f_refl - f_LOS;
                D(d).data.res(v, samp).refls(r).Delta_fD(1:no_f) = dopplerDiff; % multi-freq!
                
                %Decide material
                material = UCMP_getMaterial(currR.refls(r).objectType);
                D(d).data.res(v, samp).refls(r).material = material;
                
                Gamma = UCMP_computeReflCoeff(const, material, f_toBuild, ...
                    currR.refls(r).reflAngle_rad, rej_LHCP, currR.satPos_AER(2, prn));
                
                % Power loss
                powerLoss_ratio = abs(Gamma.Eff);
                powerLoss_dB = -20*log10(powerLoss_ratio);
                D(d).data.res(v, samp).refls(r).Delta_P(1:no_f) = powerLoss_dB;
                
                % Phase difference
                phaseDiff_refl = angle(Gamma.Eff);
                D(d).data.res(v, samp).refls(r).Delta_psi(1:no_f) = phaseDiff_refl; % multi-freq!
                
                % res(constRep, samp).fD_LOS(1:no_f,currPRN) = f_LOS - fc;
                % Doppler shift in general, needs rework
                
            end
        end
    end
end
disp('Physics computed!')

%%
% Compute distributions here, in post-processing
paramFields = {'Delta_d', 'Delta_fD', 'Delta_P', 'Delta_psi'}; % must be adapted for 3-freq
paramDenot = {'d', '\Delta f_D', '\Delta P', '\Delta \psi'};
paramNames =  {'Relative delay', 'Doppler difference', 'Power loss',  'Phase shift'};
paramUnits =   {'m', 'Hz', 'dB',  '$^\circ$'};
paramMult = [1, 1, 1,  180/pi];
for d = 1:1:length(D)
    clear paramDist
 
    for field = 1:1:length(paramFields)
        paramDist.(paramFields{field}) = zeros(3,0);
    end
    totalAdd = 0;
    for v = 1:1:size(D(d).data.res,1)
        for samp = 1:1:size(D(d).data.res,2)
            if ~isempty(D(d).data.res(v, samp).refls)
                for r = 1:1:length(D(d).data.res(v, samp).refls)
                    totalAdd = totalAdd + 1;
                    paramDist.LOS(totalAdd) =  D(d).data.res(v, samp).satLOS(D(d).data.res(v, samp).refls(r).prn);
                    paramDist.planeType{totalAdd} = D(d).data.res(v, samp).refls(r).planeType;
                    paramDist.objectType{totalAdd} = D(d).data.res(v, samp).refls(r).objectType;
                    paramDist.elev(totalAdd) = ...
                        D(d).data.res(v, samp).satPos_AER(2,D(d).data.res(v, samp).refls(r).prn);
                    paramDist.azim(totalAdd) = ...
                        D(d).data.res(v, samp).satPos_AER(1,D(d).data.res(v, samp).refls(r).prn);
                    for field = 1:1:length(paramFields)
                        newData = D(d).data.res(v, samp).refls(r).(paramFields{field});
                        paramDist.(paramFields{field}) =...
                            [paramDist.(paramFields{field}), newData'];
                    end
                end
            end
        end
    end
    D(d).paramDist = paramDist;
end


% We need for distributions:
% Total number of:
% - SPLOS
% - MP
% - NLOS
% - Blocked by building
% - Too low elevation

% MP and NLOS are combined, since a satellite can be both MP wall and MP
% ground

for d = 1:1:length(D)
    clear recModes received
    for v = 1:1:D(d).data.noV
        received(v).prns = zeros(1,0) ; %List of PRNS that are in categories 1-5 for entire simulation (for color assignment)
    end
    moreThanOneReflPerPRN = 0;
    res = D(d).data.res;
    
    D(d).totalRefls = 0;
    % Check that the reception modes have no common PRNs
    for v = 1:1:size(res,1)
        
        for samp = 1:1:size(res,2)
            %    find(strcmp({res(constRep, samp).refls.objectType } , 'build') & strcmp({res(constRep, samp).refls.planeType } , 'wall'))
            if length(res(v, samp).refls) > 0
                newNoRefls = length(find(strcmp({res(v, samp).refls.objectType } , 'build') & strcmp({res(v, samp).refls.planeType } , 'wall')));
            else
                newNoRefls = 0;
            end
            D(d).totalRefls  = D(d).totalRefls + newNoRefls;
            tooLow = find(res(v, samp).satPos_AER(2,:) < 0);
            LOS = find(res(v, samp).satLOS);
            NLOS = find(~res(v, samp).satLOS);
            
            received(v).prns = union(received(v).prns, LOS);
            
            anyRefls = [];
            refls_per_PRN = zeros(1,31);
            for r = 1:1:res(v, samp).noRefl
                if strcmp(res(v,samp).refls(r).objectType,'build') && strcmp(res(v,samp).refls(r).planeType, 'wall')
                    currPRN = res(v,samp).refls(r).prn;
                    anyRefls = [anyRefls, currPRN];
                    refls_per_PRN(currPRN) = refls_per_PRN(currPRN) + 1;
                end
            end
            
            if any(refls_per_PRN > 1)
                %   fisk
                moreThanOneReflPerPRN = moreThanOneReflPerPRN + 1;
            end
            received(v).prns = union(received(v).prns, anyRefls);
            
            recModes(v, samp).SPLOS = setdiff(LOS, anyRefls);
            recModes(v, samp).MP = intersect(anyRefls, LOS);
            recModes(v, samp).NLOS = intersect(anyRefls, NLOS);
            recModes(v, samp).Blocked = setdiff(setdiff(NLOS, tooLow), anyRefls);
            recModes(v, samp).tooLow = tooLow;
            
        end
    end
    
    if moreThanOneReflPerPRN > 0
        disp('Careful: Some satellites were reflected via several paths')
    end
    
    D(d).received = received;
    D(d).recModes = recModes;
end
disp('Distributions computed!')
%% Visualiation of reflections
v = 1 % Which of the vehicles should be plotted?
fontSize = 40
simToPlot = 1
res = D(simToPlot).data.res
P = D(simToPlot).data.P
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
for w = 1:1:length(P)-5
    poly_rectangle(P(w).points(:,1), P(w).points(:,2),...
        P(w).points(:,3), P(w).points(:,4),  blockColors.(P(w).objectType))
    
    % Add plane ID number, for debuggung purposes (optional)
    %     text(P(w).lims(1,1) + (P(w).lims(1,2) -  P(w).lims(1,1))/2, ...
    %         P(w).lims(2,1) + (P(w).lims(2,2) -  P(w).lims(2,1))/2, ...
    %         P(w).lims(3,1) + (P(w).lims(3,2) -  P(w).lims(3,1))/2, ...
    %         num2str(w), 'FontSize', 15)
end

view(0, 0)
%view(-10, 43)
xlabel('East (m)' , 'FontSize', fontSize, 'Interpreter', 'latex')
ylabel('North (m)', 'FontSize', fontSize, 'Interpreter', 'latex')
zlabel('Up (m)', 'FontSize', fontSize, 'Interpreter', 'latex')

zoomPoint = [0 0 0]

% Option to zoom on car initial position:
zoomSize = 250; % m
%centerpoint = a.pos_init
% xlim([zoomPoint(1)-zoomSize/2    zoomPoint(1)+zoomSize/2])
% ylim([zoomPoint(2)-zoomSize/2    zoomPoint(2)+zoomSize/2])
%zlim([0 100])

disp('Urban canyon plotted!')
xlim([0 D(simToPlot).data.upperPlotLim])
ylim([0 D(simToPlot).data.upperPlotLim])
% zlim([0 D(simToPlot).data.upperPlotLim])

xlim([-150 150])
ylim([-150 150])
zlim([0 250])


ax = gca;
ax.FontSize = fontSize;
ax.TickLabelInterpreter = 'latex';

title(['Reflections for collaborative positioning'])
%%
% Plot reflection(s)
disp('Plotting signal paths...')
no_PRNs = length(D(d).received(v).prns);
% Received PRNS must depebd on startIt!
colsToUse = possibleCols(1:1:no_PRNs);
hold on
for samp = 1:1:size(res,2)
    % Plot satellite position, antenna position, LOS
    scatter3(res(v, samp).antPos_ENU(1), res(v, samp).antPos_ENU(2),...
        res(v, samp).antPos_ENU(3), 100, 'k', 'Marker', 'o') % Vehicle position
    
    %     for p = s.allPRN
    %         scatter3(res(step).spos(1,p), res(step).spos(2,p), res(step).spos(3,p), 1000, 'k') % Satellite position
    %     end
    
    for r = 1:1:length(D(simToPlot).recModes(samp).SPLOS)
        if length(D(simToPlot).recModes(samp).SPLOS) == 0
            break
        end
        currPRN = D(simToPlot).recModes(samp).SPLOS(r);
        colInd =  find(currPRN==received(v).prns);
        currCol = colsToUse{colInd};
        currLOSvec =  [res(v, samp).satPos_ENU(:,currPRN),  res(v, samp).satArrival(:,currPRN)]'; % maybe switch
        
        plot3(currLOSvec(:,1), currLOSvec(:,2), currLOSvec(:,3),...
            'color', currCol, 'LineStyle', '-') % Line of sight (can also be blocked)
    end
    
    for r = 1:1:res(v, samp).noRefl %%% Here must not only reflections be plotted, also SPLOS!!!
        
        currPRN = res(v, samp).refls(r).prn;
        colInd =  find(currPRN==received(v).prns);
        currCol = colsToUse{colInd};
        currLOSvec = [res(v, samp).satPos_ENU(:,currPRN),  res(v, samp).satArrival(:,currPRN)]'; % maybe switch
        
        if res(v, samp).satLOS(currPRN) == 1
            plot3(currLOSvec(:,1), currLOSvec(:,2), currLOSvec(:,3),...
                'color', currCol, 'LineStyle', '-') % Line of sight (can also be blocked)
        end
        reflPoint = res(v, samp).refls(r).reflPoint;
        %mirr = res(step).refls(r).mirrPoint;
        %proj = res(step).refls(r).projPoint;
        
        %         D = v.pos - s.pos;
        %         B = refl - s.pos;
        Bplot = [reflPoint,  res(v, samp).satPos_ENU(:,currPRN)]';
        Rplot = [res(v, samp).antPos_ENU, reflPoint]';
        
        hold on
        scatter3(reflPoint(1), reflPoint(2), reflPoint(3), 'MarkerFaceColor', currCol, 'Marker', '*')
        %scatter3(mirr(1), mirr(2), mirr(3), 'b')
        %scatter3(proj(1), proj(2), proj(3), 'y')
        
        plot3(Bplot(:,1), Bplot(:,2), Bplot(:,3), 'color', currCol, 'LineStyle', '-.')
        plot3(Rplot(:,1), Rplot(:,2), Rplot(:,3), 'color', currCol, 'LineStyle', '-.')
        
    end
end

clear prnLeg
h = zeros(1, length(received(v).prns));
for r = 1:1:length(received(v).prns)
    h(r) = plot3(NaN,NaN,NaN, 'color', colsToUse{r});
    prnLeg{r} = ['PRN ' num2str(received(v).prns(r))];
end
legend(h, prnLeg , 'FontSize', fontSize, 'Interpreter', 'latex');

disp('Signal paths plotted!')
%%
% Simulate receivers to compute pseudorange errors and resulting SNR (2 frequencies)

% Receiver definition
sampFactor = 0.1 %24 % sampFactor/12 must be intenger for E5!
corr_spac = 0.5 % spacing in EPL correlator, number of chips
noChips_del = 2 % Number of chips in each direction to find correlation peak
noChips_SNR = 30 % Number of chips in each direction to decide noise level for SNR
noise_phase = 0 % phase noise of receiver
noise_amp = 0 % amplitude noise of receiver
cb = 1e3; % clock bias of receiver (not yet used)
sys = 'GPS'
freqBands = {'L1', 'L2'}
comps = {'CA', 'CM'}
no_freqs = 1:1:length(freqBands);

% Compute pseudorange errors and resultin SNR
clear obs_with_mp
lines  = 0;
disp('Starting receiver simultion to compute pseudorange errors and new SNRs...')
tic
for v = 1:1:size(res,1)
    for s = 1:1:size(res,2)
        s
        clear okInds
        for r = 1:1:length(res(v,s).refls) % exclude reflections on car
            okInds(r) = ~strcmp(res(v,s).refls(r).objectType, 'car');
        end
        currRefls = res(v,s).refls(okInds);
        allPRNs = [currRefls.prn];
        [C,IA,IC] =  unique(allPRNs);
        
        for sat = 1:1:length(IA) % Loop over all satellites with at least 1 reflection
            % Here we create correlation for one vehicle, one time
            % step, and all received signals from a given PRN
            
            satID = C(sat);
            newPRNinds = IA(sat):1:(IA(sat) + sum(IC==sat) - 1); % Store all reflections here

            % Check if  LOS signal is received
            LOS_received = res(v,s).satLOS(satID);

            noRefls = length(newPRNinds);
            clear delta_d delta_P
            for refl = 1:1:noRefls % Loop over invidiual reflections
                delta_d(refl,1:2) = currRefls(newPRNinds(refl)).Delta_d(1:no_freqs);
                delta_P(refl,1:2) = currRefls(newPRNinds(refl)).Delta_P(1:no_freqs);
            end
            
            [prErr, SNR_factor] = compPRerr(sampFactor, corr_spac, noChips_del, ...
                noChips_SNR, noise_phase, noise_amp, cb, sys, freqBands, comps, satID,...
                LOS_received, noRefls, delta_d, 10.^(-delta_P/20), const);
            
            prDiff_L1 = prErr(1);
            prDiff_L2 = prErr(2);
            
            % Elevation-based SNR model: https://www.sciencedirect.com/science/article/pii/S1110982316300412
            elev = res(v,s).satPos_AER(2,satID);
            SNR_L1 = 3.199e-5 * elev ^3 -0.0081*elev ^2 +0.6613*elev + 31.38;
            SNR_L2 = SNR_L1 - 3;
            
            % Computed new SNR according to receiver output:
            SNR_MP_L1 = 20*log10( 10.^(SNR_L1/20) * SNR_factor(1));
            SNR_MP_L2 = 20*log10( 10.^(SNR_L2/20) * SNR_factor(2));
            
            lines = lines + 1;
            % Save results (to be exported to obs files)
            obs_with_mp(v,lines).data = [...
                res(v,s).currTime, ... % Time
                satID,... % Satellite PRN
                prDiff_L1,...  % This is to be added to frequency 1 pseudorange
                SNR_MP_L1, ... % This is the new SNR on frequency 1
                prDiff_L2, ... % This is to be added to frequency 2 pseudorange
                SNR_MP_L2... % This is the new SNR on frequency 2
                ];

        end

    end
end
receiverSimTime = toc
disp('New mutlipath-influenced pseudoranges and SNRs computed!')
% Missing: Here the Rinex obs file have to be modified with new
% pseudoranges and SNR values for each satellite and frequency
% Unfinished code for Rinex export in functions/exportRinexObs.m