function [refl, mirror, proj] = UCMP_computeReflection(normal, planePoint, planeLims, ant, sat)
%COMPUTEREFLECTION Summary of this function goes here
%   Detailed explanation goes here
         
% normal = P(w1).normal
% planePoint = P(w1).points(:,1)
% planeLims =  P(w1).lims
% ant = antPos_ENU
% sat =satPos_ENU

proj = ant + sum((planePoint-ant ).* normal) / sum(normal .* normal) * normal;
mirror = ant + 2*(proj-ant);

refl = UCMP_planeLineIntersect(normal, planePoint, planeLims, mirror, sat);

end

