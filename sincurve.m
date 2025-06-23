function [ysin] = sincurve(ymin,ymax,ndays,tmax)

%........................................................................
T = ndays;
tau = ndays/4;
time = [1:1:tmax];
ttime = time - tau;
wrad = (2*pi) / T;
%........................................................................
A0 = ymin;
A = (ymax-ymin)/2; %amplitud.
ysin = A0 + A*( sin(wrad*ttime) + 1); %Sinusoidal function.
%........................................................................


return