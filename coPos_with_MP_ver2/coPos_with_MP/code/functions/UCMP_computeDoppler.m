function f_out = UCMP_computeDoppler(p_rec, v_rec, p_emi, v_emi, f_in, const)
%COMPUTEDOPPLER Summary of this function goes here
%   Detailed explanation goes here
% Computes fD as difference between carrier frequency and received
% frequency

% https://isaacphysics.org/concepts/cp_doppler_effect
%    https://en.wikipedia.org/wiki/File:DopplerSatScheme.png
deltav = v_emi - v_rec;
carsatvec = p_emi - p_rec;
% doppAngle = acos(       dot(deltav, carsatvec) / (norm(deltav) * norm(carsatvec))   );
% f_out = f_in * ( 1 + norm(deltav)/const.c *cos(doppAngle));

 f_out = f_in * ( 1 + sum(deltav.*carsatvec)/ ( norm(carsatvec) * const.c   ) );

end

