% [obscircle1] = readRinexObsSprient("/home/liy4hi/000/MASTERTHESIS/COOOOOOODE/meas/Obs/circle1.obs")
% [obscircle2] = readRinexObsSprient("/home/liy4hi/000/MASTERTHESIS/COOOOOOODE/meas/Obs/circle2.obs")
[obscircle] = readRinexObsSprient("/home/liy4hi/000/MASTERTHESIS/COOOOOOODE/meas/Obs/circle.obs");
 car = CarModel("/home/liy4hi/000/MASTERTHESIS/COOOOOOODE/meas/Car/circle.umt", obscircle);
% car2 = CarModel("/home/liy4hi/000/MASTERTHESIS/COOOOOOODE/meas/Car/circle2.umt", obscircle2);
car1 = load("/home/liy4hi/000/MASTERTHESIS/COOOOOOODE/meas/Car/circle1.umt");
car2 = load("/home/liy4hi/000/MASTERTHESIS/COOOOOOODE/meas/Car/circle2.umt");
r = importdata("values.csv")

%%
[startpoint1(1),startpoint1(2),startpoint1(3)] = ecef2geodetic(wgs84Ellipsoid('meters'),car1(1,3),car1(1,4),car1(1,5));

[car1ENU(:,1),car1ENU(:,2),car1ENU(:,3)] = ecef2enu(car1(:,3),car1(:,4),car1(:,5),...
                              startpoint1(1),startpoint1(2),startpoint1(3),wgs84Ellipsoid);

[startpoint2(1),startpoint2(2),startpoint2(3)] = ecef2geodetic(wgs84Ellipsoid('meters'),car1(1,3),car1(1,4),car1(1,5));

[car2ENU(:,1),car2ENU(:,2),car2ENU(:,3)] = ecef2enu(car2(:,3),car2(:,4),car2(:,5),...
                              startpoint2(1),startpoint2(2),startpoint2(3),wgs84Ellipsoid);

figure(1)
plot(car1ENU(:,1),car1ENU(:,2),'.')
hold on 
grid on
plot(car2ENU(:,1),car2ENU(:,2),'*')
axis equal
legend('Vehicle1','Vehicle2')
xlabel('East West direction (meter)')
ylabel('North South direction (meter)')

[rENU(:,1),rENU(:,2),rENU(:,3)] = ecef2enu(r(:,1),r(:,2),r(:,3),...
                              startpoint1(1),startpoint1(2),startpoint1(3),wgs84Ellipsoid);

rENU(1801,:) = [];

figure(2)
plot(car1ENU(:,1)-rENU(:,1), 'r');
hold on
plot(car1ENU(:,2)-rENU(:,2), 'b');
hold on
grid on
plot(car1ENU(:,3)-rENU(:,3), 'g');
legend('X','Y','Z');
xlabel('Time stamp');
ylabel('Position difference (meter)');



