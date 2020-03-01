function [ output_args ] = UCMP_generateCorrOutput( input_args )
%UCMP_GENERATECORROUTPUT Summary of this function goes here
%   Detailed explanation goes here

init ; 
sys = 'GPS'
freqBand = 'L1'

satID = 1

f_samp = const.GPS.L1.CA.f_carrier*20;
f_chip = const.GPS.L1.CA.f_chip
f_carr =  const.GPS.L1.CA.f_carrier;

codeLength_s = 1e-3

no_samps = f_carr/f_chip % 1540 complete carrier waves fit inside the complete code

tVec = 0:codeLength_s/(no_samps-1):codeLength_s;


codes_upsampled = codeGenerator(sys, freqBand, satID, f_samp, const);
size(codes_upsampled.GPS.L1.CA)

LOSsig = cos(2*pi*f_carr*tVec); % * code


DeltaP = 3 % dB
DeltaD =  20 % m
DeltafD = 0 % Hz
DeltaPhi = pi %

refl = 10^(-DeltaP/20) * cos( 2*pi*(f_carr + DeltafD)*(tVec + DeltaD/const.c) + DeltaPhi) % * shifted code


figure;
plot(tVec, LOSsig, 'g')
hold on
plot(tVec, refl, 'r')
plot(tVec, refl+LOSsig, 'b')



end

