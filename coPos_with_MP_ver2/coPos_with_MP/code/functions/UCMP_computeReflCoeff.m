function Gamma = UCMP_computeReflCoeff(const, material, freqList, reflAngle, rej_LHCP, elevAngle) %% ADD ELEV???
%COMPUTE_ Summary of this function goes here
%   Detailed explanation goes here
%

epsil_r = const.relPerm.(material);
sigma = const.cond.(material);
lambdas = const.c ./ freqList; % Wavelengths

% epsilon: relative permittivity
% sigma: dielectric constant 

epsil = epsil_r - 1j*60*lambdas*sigma; % Complex dielectric constants
Gamma.Horizontal = (sin(reflAngle) - sqrt(epsil - cos(reflAngle).^2)) ./ (sin(reflAngle) + sqrt(epsil - cos(reflAngle).^2)); % Horizontal reflection coefficient
Gamma.Vertical = (epsil*sin(reflAngle) - sqrt(epsil - cos(reflAngle).^2)) ./ (epsil*sin(reflAngle) + sqrt(epsil - cos(reflAngle).^2)); % Vertical reflection coefficient
Gamma.Copolar = (Gamma.Horizontal + Gamma.Vertical) / 2;
Gamma.Crosspolar = (Gamma.Horizontal - Gamma.Vertical) / 2 ;
Gamma.Eff = (abs(Gamma.Copolar) + 10^(-rej_LHCP*sind(elevAngle)/20)*abs(Gamma.Crosspolar))*exp(-j*pi); % Effective reflection coefficient

end

