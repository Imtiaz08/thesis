function [S_carrier, theta_MP] = genGALE1(code, f_samp, t, amp, std_noise_phase, std_noise_amp, time_delay_MP, const)
%GENGALE1 Summary of this function goes here
%   Detailed explanation goes here
% Constants
%Dperiod = 20e-3;
f_carrier = const.GAL.E1.f_carrier;
f_subcarrier = const.GAL.E1.f_subcarrier;

%noSamples = length(t);
theta_MP_mc = 2*pi*mod(f_carrier*time_delay_MP, 1);
samples_delay_MP = round(time_delay_MP*f_samp);

code_E1B = circshift(code(1,:), samples_delay_MP);
code_E1C = circshift(code(2,:), samples_delay_MP);

data_E1B = 1;

% baseband of I component
alpha = sqrt(10/11);
beta = sqrt(1/11);
m = 6; % 6
n = 1; % 1
theta_MP_sc1 = 2*pi*mod(1*f_subcarrier*time_delay_MP, 1);
theta_MP_sc6 = 2*pi*mod(6*f_subcarrier*time_delay_MP, 1);

sc_E1_Ba = sign(sin(2*pi* n*f_subcarrier*t+  std_noise_phase*randn(1,length(t)) + theta_MP_sc1 ));
sc_E1_Bb = sign(sin(2*pi* m*f_subcarrier*t+ std_noise_phase*randn(1,length(t)) + theta_MP_sc6));
% sc_Ei_Ca = sign(sin(2*pi* n*f_chip*t+ phi*randn(1,length(t))+theta_MP_sc1 ));
% sc_Ei_Cb = sign(sin(2*pi* m*f_chip*t+ +phi*randn(1,length(t))+theta_MP_sc6));

Bcomp = (code_E1B.*data_E1B).*(alpha*sc_E1_Ba + beta*sc_E1_Bb);
Ccomp = code_E1C.*1.*(alpha*sc_E1_Ba - beta*sc_E1_Bb);

% tToPlot =  1/f_chip; Verify plots from definition
% indToPlot = 1:1:tToPlot*f_samp;
% figure;
% plot(t(indToPlot), Bcomp(indToPlot),'b')
% figure;
% plot(t(indToPlot), Ccomp(indToPlot),'r')

I_bb = 1/sqrt(2).*(Bcomp - Ccomp);
%Q_bb = zeros(1, length(I_bb)); % E1 A component empty
%I_carrier = .*E1_total;
%s_bb = I_bb + 1j*Q_bb;

S_carrier = sqrt(2*amp)*(I_bb.*cos(2*pi*f_carrier*t + std_noise_phase*randn(1,length(t)) + theta_MP_mc)) +...
    std_noise_amp*randn(1, length(t));

theta_MP = [theta_MP_mc, theta_MP_sc1, theta_MP_sc6];
% +Q_bb.*sin(2*pi*f_carrier*t + phi + theta_MP) );
%IQ_carrier = I_carrier + 1j*Q_carrier;
% size(Q_carrier)
% size(I_carrier)
% 
% I_bb = I_carrier .*cos(2*pi*f_GAL_E1*t);
% Q_bb = Q_carrier .*sin(2*pi*f_GAL_E1*t);
% IQ_bb = I_bb +1j*Q_bb; 

end

