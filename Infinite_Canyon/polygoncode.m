function polycode = polygoncode(plane_code, polys)

polycodes_all = [11, 10, 21, 20, 31, 30, 41, 40];


for i = 1:8    
    if polys(:,:,i) == plane_code
        polycode = polycodes_all(i);
        break
    end
end
