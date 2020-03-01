function P = UCMP_generateBlock(under, upper, blockType)
%GENERATEBLOCK Summary of this function goes here
%   Detailed explanation goes here

w = 0;
% http://www.matrixlab-examples.com/3d-polygon.html
p1 = [under(1),  under(2),    under(3)]';
p2 = [upper(1),  under(2),    under(3)]';
p3 = [upper(1),  upper(2),    under(3)]';
p4 = [under(1),  upper(2),    under(3)]';
p5 = [under(1),  under(2),    upper(3)]';
p6 = [upper(1),  under(2),    upper(3)]';
p7 = [upper(1),  upper(2),    upper(3)]';
p8 = [under(1),  upper(2),    upper(3)]';

myWalls(1).points = [p5, p6, p7, p8]; % roof
% 1-4: walls
myWalls(2).points = [p2, p6, p5, p1]; % opposite myWalls(4)
myWalls(3).points = [p2, p6, p7, p3]; % opposite myWalls(5)
myWalls(4).points = [p3, p4, p8, p7];
myWalls(5).points = [p1, p4, p8, p5];

myWalls(6).points = [p1, p2, p3, p4]; % floor

oppTable(1) = 6;
oppTable(2) = 4;
oppTable(3) = 5;
oppTable(4) = 2;
oppTable(5) = 3;
oppTable(6) = 1;

for m = 1:1:length(myWalls)
    
    wp1 = myWalls(m).points(:,1);
    wp2 = myWalls(m).points(:,2);
    wp3 = myWalls(m).points(:,3);
    wp4 = myWalls(m).points(:,4);
    normal = cross(wp1-wp2, wp1-wp3);
    normal = normal/norm(normal);
    
    % Force normal outwards:
    oppPoint = myWalls(oppTable(m)).points(:,1);
    % Vector from opposite wall to current wall, should align with normal:
    oppVec = wp1 - oppPoint;
    
    if abs(acos( sum(oppVec.* normal) / (norm(oppVec) * norm(normal)))) > pi/2
        normal = -1*normal; % invert direction of normal
    end
    
    allwp = [wp1, wp2, wp3, wp4];
    lims = [min(allwp(1,:))   max(allwp(1,:));...
        min(allwp(2,:))   max(allwp(2,:));...
        min(allwp(3,:))   max(allwp(3,:))];
    
    % Add new planes to list: (ignore floors)
    if m < 6
        w = w + 1;
        P(w).normal = normal;
        P(w).lims = lims;
        P(w).points = [wp1, wp2, wp3, wp4];
        
        % Specifiy wall type
        if m == 1
            P(w).planeType = 'roof';
        elseif m == 6
            P(w).planeType = 'floor';
        else
            P(w).planeType = 'wall';%['wall-' num2str(m)];
        end
        P(w).objectType = blockType;
    end
end


end

