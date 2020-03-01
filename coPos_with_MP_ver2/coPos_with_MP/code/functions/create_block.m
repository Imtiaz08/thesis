function create_block(center, x_size, y_size, z_size)
%CREATE_BLOCK Summary of this function goes here
%   Detailed explanation goes here

p1 = center + [0 0 0];
p2 = center + [x_size 0 0];
p3 = center + [x_size y_size 0];
p4 = center + [0 y_size 0]; 

p5 = p1 + [0 0 z_size];
p6 = p2 + [0 0 z_size];
p7 = p3 + [0 0 z_size];
p8 = p4 + [0 0 z_size];

x = [p1(1) p2(1) p3(1) p4(1)];
y = [p1(2) p2(2) p3(2) p4(2)];
z = [p1(3) p2(3) p3(3) p4(3)];

fill3(x, y, z, 1);

x = [p5(1) p6(1) p7(1) p8(1)];
y = [p5(2) p6(2) p7(2) p8(2)];
z = [p5(3) p6(3) p7(3) p8(3)];
fill3(x, y, z, 2);

x = [p2(1) p6(1) p7(1) p3(1)];
y = [p2(2) p6(2) p7(2) p3(2)];
z = [p2(3) p6(3) p7(3) p3(3)];
fill3(x, y, z, 3);
 
x = [p2(1) p6(1) p5(1) p1(1)];
y = [p2(2) p6(2) p5(2) p1(2)];
z = [p2(3) p6(3) p5(3) p1(3)];
fill3(x, y, z, 4);

poly_rectangle(p3, p4, p8, p7)
poly_rectangle(p1, p4, p8, p5)

end

