init

fileName = 'testRin'
folder = 'U:/Minimization_GNSS_Errors/code'

path = [folder '/' fileName '.obs']


% > 2019 07 03 12 00  0.0000000  0  8       0.000000000000
% G07  25371998.322 9 133330704.775 9      2213.624 9        54.123 9  25371999.005 6         0.000 6      1724.902 6        37.977 6
% G10  23849385.600 9 125329311.920 9      3537.300 9        54.655 9  23849389.366 6         0.000 6      2756.338 6        37.977 6

year_start = 2019
month_start = 7
day_start = 3

hour_start = 12
minute_start = 0
second_start = 1

vecId = 1
samp = 1

% General parameters for this vehicle and sample
rec(vecId).meas(samp).time = 100
rec(vecId).meas(samp).no_obs = 2
rec(vecId).meas(samp).unknown_1 = 0.000000000000
% Per satellite counter
satCounter = 1
rec(vecId).meas(samp).sats(satCounter).satId = 7
rec(vecId).meas(samp).sats(satCounter).sys = 'GPS'
% freq 1
rec(vecId).meas(samp).sats(satCounter).rho_f1 = 25371998.322
rec(vecId).meas(samp).sats(satCounter).unknown_1 = 9
rec(vecId).meas(samp).sats(satCounter).cp_f1 = 133330704.775
rec(vecId).meas(samp).sats(satCounter).unknown_2 = 9
rec(vecId).meas(samp).sats(satCounter).fD = 2213.624
rec(vecId).meas(samp).sats(satCounter).unknown_3 = 9
rec(vecId).meas(samp).sats(satCounter).CN0 = 54.123
rec(vecId).meas(samp).sats(satCounter).unknown_4 = 9
% freq 2
rec(vecId).meas(samp).sats(satCounter).rho_f2 = 25371999.005
rec(vecId).meas(samp).sats(satCounter).unknown_1 = 6
rec(vecId).meas(samp).sats(satCounter).cp_f1 = 0.000
rec(vecId).meas(samp).sats(satCounter).unknown_2 = 6
rec(vecId).meas(samp).sats(satCounter).fD = 1724.902
rec(vecId).meas(samp).sats(satCounter).unknown_3 =6
rec(vecId).meas(samp).sats(satCounter).CN0 = 37.977
rec(vecId).meas(samp).sats(satCounter).unknown_4 = 6

% > 2019 07 03 12 00  0.0000000  0  8       0.000000000000
timeHeader_format = {'%01.f', '%02.f', '%02.f', '%02.f', '%02.f', '%.7f',...
    '%01.f' '%01.f''%.12f' }
timeHeader_spaces = [' ', ' ', ' ']
% var name, format

% types:
% integer with 1 leading 0
% integer with no leading 0
% double with 3 decimals
% double with 7 decimals
% double with 12 decimals

test1 = 7
test2 = 11
num2str(test2,'%02.f')
num2str(test1, '%01.f')
datastring=num2str(1231232.129830184598230523543,'%.3f')
datastring=num2str(2.129830184598230523543,'%.7f')
datastring=num2str(2.129830184598230523543,'%.12f')