function [D] = diffprofile(Ddeep,Dsurface,mld,zmin,zmax,dz)
% after Roger Cropp Griffith University, Brisbane; 2004 DMS modeling

% abs_deltaZ = 10; % Weak slope
abs_deltaZ = 20; % Midium slope (USE THIS ONE!!!)
% abs_deltaZ = 40; % Strong slope

if Dsurface > Ddeep % For "D" being turbulent diffusion or water temperature (which decrease with depth)
    Dmax = Dsurface;
    Dmin = Ddeep;
    deltaZ = abs_deltaZ * (-1); 
else %% elseif Dsurface <= Ddeep % For "D" being water density (which increases with depth)
    Dmax = Ddeep;
    Dmin = Dsurface;
    deltaZ = abs_deltaZ;
end

% SIGMOIDAL EQUATION FOR DIFFUSION PROFILE:
zaxis = [zmin:dz:zmax-dz];
Dstep = (deltaZ / (Dmax - Dmin));
A = Dstep * (Dmin - Dmax) * 2 * (zaxis - mld) / zmax;
D = (Dmax + Dmin .* exp(A)) ./ (1 + exp(A));


return