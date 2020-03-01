function [refl, mirrpoint, polycode, flag] = reflectionline(L1_point,L1_poly1, L1_poly2, L2_point, L2_interpoly, L2_noninter, L2_decide, otherpoly, sat_pos, cross, L_cross, polys)
% L1_point = symPO_n;
% L1_poly1 = poly1_E;
% L1_poly2 = poly2_E;
% 
% L2_point = symPO_e;
% L2_interpoly = poly4_N;
% L2_noninter = poly3_N;
% L2_decide = poly4_E;
% otherpoly = poly3_E;
% L1_point is possible mirror point1, L2_point is possible mirror point2

refl = 5;
flag = [];
reflectionpoint = zeros(1,3);
mirrpoint = zeros(1,3);
polycode = 000;

L1 = createLine3d(L1_point, sat_pos);                        
PO1_L1 = intersectLinePolygon3d(L1, L1_poly1);
PO2_L1 = intersectLinePolygon3d(L1, L1_poly2);

L2 = createLine3d(L2_point, sat_pos);
PO_interpolyL2 = intersectLinePolygon3d(L2, L2_interpoly);
PO_noninterL2 = intersectLinePolygon3d(L2, L2_noninter);

if 1 && ismember(0,isnan([PO1_L1, PO2_L1]))   % intersect with plane N               
    PO3_E = intersectLinePolygon3d(L1, otherpoly);
    PO3_N = intersectLinePolygon3d(L1, L2_noninter);
    PO4_E = intersectLinePolygon3d(L1, L2_decide);
    PO4_N = intersectLinePolygon3d(L1, L2_interpoly);

    if ismember(1, isnan([PO3_E, PO3_N, PO4_E, PO4_N])) == 0   % intersect with other polys
        refl = 0;
    else
        refl = 1; 
        flag = [flag; 1];
        if 1 && all(isnan(PO1_L1)==1)  % no point with PO1_E
            polycode = polygoncode(L1_poly2, polys);
            reflectionpoint = PO2_L1; 
            mirrpoint = L1_point;
        else
            polycode = polygoncode(L1_poly1, polys);
            reflectionpoint = PO1_L1;
            mirrpoint = L1_point;
        end
    end
end 


if 1 && all(isnan(PO_interpolyL2)==0) && all(isnan(PO_noninterL2) == 1)       % intersect with L2_interpoly not poly 3_N
    PO_decideL2 = intersectLinePolygon3d(L2, L2_decide);
    if abs(PO_decideL2(1) - cross(1)) > L_cross
        refl = 1;
        polycode = polygoncode(L2_interpoly, polys);
        reflectionpoint = PO_interpolyL2;
        mirrpoint = L2_point;
        flag = [flag; 1];
    else
        refl = 0;
    end
end
