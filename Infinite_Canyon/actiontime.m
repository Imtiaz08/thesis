function [c] = actiontime(time, t0)

% get time in format of 00:00:00 

delta = time-t0;

hour = floor(delta/3600);
mint = floor((delta-hour*3600)/60);
sec = mod(delta, 60);

T = [hour, mint, sec];
D = duration(T);
c = char(D,'hh:mm:ss');
