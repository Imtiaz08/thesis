function crossPoint = UCMP_planeLineIntersect(normal, planePoint, planeLims, linePoint1, linePoint2)
%PLANELINEINTERSECT Summary of this function goes here
%   Detailed explanation goes here
% normal = P(w).normal
%   planePoint           =   P(w).points(:,1)
% planeLims = P(w).lims
% 
% 
% linePoint1 = antPos_ENU
% linePoint2 = satPos_ENU

% Intersect plane and line
crossPoint = linePoint2 + sum((planePoint-linePoint2).* normal) /...
    sum((linePoint1 - linePoint2).* normal) * (linePoint1 - linePoint2); 

sens = 1e-3; % 1 mm
% Check that crosssPoint is within plane limits
crossPointWithinPlane = NaN*ones(3,1);
for d = 1:1:3
    crossPointWithinPlane(d) = crossPoint(d) >= planeLims(d,1) - sens &&...
        crossPoint(d) <= planeLims(d,2) + sens;
    if ~ crossPointWithinPlane(d) % Optimiziation
        crossPoint = NaN*ones(3,1);
        return
    end
end

% if ~all(crossPointWithinPlane)
%     crossPoint = NaN*ones(1,3);
%     return
% end


% Check that cross point lies between satellite and antenna
v = (linePoint2 - linePoint1);
t = (crossPoint - linePoint1)' / v';
crossIsPointBetweenAntSat = (t >= 0 && t <= 1);

if ~crossIsPointBetweenAntSat
    crossPoint = NaN*ones(3,1);
end

% linePoints = [linePoint1; linePoint2]
% underLine = min(linePoints)
% upperLine = max(linePoints)
% betweenX = crossPoint(1) >= underLine(1) && crossPoint(1) <= upperLine(1);
% betweenY = crossPoint(2) >= underLine(2) && crossPoint(2) <= upperLine(2);
% betweenZ = crossPoint(3) >= underLine(3) && crossPoint(3) <= upperLine(3);
% 
% 
% 
% if ~(withinX && withinY && withinZ && betweenX && betweenY && betweenZ)
%     crossPoint = NaN*ones(1,3);
% end

end

