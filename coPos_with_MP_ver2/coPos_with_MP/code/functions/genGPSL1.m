function [S_carrier, theta_MP] = genGPSL1(code, f_samp, t, Pc, std_noise_phase, std_noise_amp, time_delay_MP, const)
%GENGPSL1 Summary of this function goes here
%   Detailed explanation goes here
 
% Constants
% CAperiod = 1e-3;
% Dperiod = 20e-3;

f_carrier = const.GPS.L1.f_carrier;
%XPperiod = 1/10.23e6;

%Pp = Pc/(10^(3/10));

theta_MP = 2*pi*mod(f_carrier*time_delay_MP, 1);
samples_delay_MP = round(time_delay_MP*f_samp);

%noSamples = length(t);


% Generate navigation data message (fake random)
%noDbits = stopTime/Dperiod; % Total number of data bits to fit in time frame
%Dbits = randi([0 1],1, ceil(noDbits)); % Generate base data bits
% Dbits = 1*ones(1, ceil(noDbits));
% Dbits(Dbits==0) = -1;
% D = repmat(Dbits,[1, ceil(noSamples/length(Dbits))]);
D = 1;

% Generate Gold Code
%XG = repmat(code,[1 length(t)/length(CAcode)]);
XG = circshift(code, samples_delay_MP);

% Generate precision code (fake random)
% noXPbits = stopTime/XPperiod;
% XPbits = randi([0 1],1, ceil(noXPbits ));
% XPbits(XPbits==0) = -1;
% XP = repmat(XPbits,[1, ceil(noSamples/length(XPbits))]);

I_bb = sqrt(2*Pc) .*XG .* D;
%Q_bb = zeros(1,length(I_bb)); %sqrt(2*Pp) .* XP .* D; simulate Q to empty
%s_bb = I_bb + 1j*Q_bb;
% Time version
I_carrier = 1/sqrt(2) * I_bb .*cos(2*pi*f_carrier*t + std_noise_phase*randn(1,length(t)) + theta_MP);
% Q_carrier = Q_bb .* sin(2*pi*f_carrier*t +phi*randn(1,length(t)) + theta_MP);
S_carrier = I_carrier + std_noise_amp*randn(1, length(t));

% Baseband version
% I_bb = I_carrier .*cos(2*pi*f_GPS_L1*t);
% Q_bb = Q_carrier .*sin(2*pi*f_GPS_L1*t);
% IQ_bb = I_bb +1j*Q_bb; 
end

