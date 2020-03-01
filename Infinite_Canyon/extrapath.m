function [ distance ] = extrapath(car, sat, mirrorpoint)

line = createLine3d(mirrorpoint, sat);
pjPO = projPointOnLine3d(car, line);
distance = norm([mirrorpoint - pjPO]);