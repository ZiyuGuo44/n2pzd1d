function [DIFF] = diffusion_i(Cz,dt,dz,kz,nz,keyScheme)

%[ see "SPEAD_1D_TurbulentDiffusion.m" coded by Guillaume Le Gland ]
%Add implicit mode to avoid instability, as in NEMO 
%When time step is much larger than dz^2/(2*KZ), the explicit scheme is
%unstable, the Crank-Nicolson scheme creates oscillations and the implicit
%scheme homogenizes the concentrations. This is why we need implicit.
% For now, only works with Neumann boundary condition

T = zeros(nz); % Inverse of the transport matrix

if strcmp(keyScheme,'Implicit')

    for i = 2:nz-1
        T(i-1,i) = - kz(i)   * dt / dz^2;
        T(i+1,i) = - kz(i+1) * dt / dz^2;
        T(i,i) = 1.0 - T(i-1,i) - T(i+1,i);
    end
    T(2,1)     = - kz(2)  * dt / dz^2;
    T(1,1)     = 1.0 - T(2,1);
    T(nz-1,nz) = - kz(nz) * dt / dz^2;
    T(nz,nz)   = 1.0 - T(nz-1,nz);

    % Implicit time derivative
    DIFF = (T \ Cz - Cz) / dt;

elseif strcmp(keyScheme,'Explicit')
    
    % split = 500;
    split = 100;
    for i = 2:nz-1
        T(i-1,i) = kz(i)   * dt / (split*dz^2);
        T(i+1,i) = kz(i+1) * dt / (split*dz^2);
        T(i,i) = 1.0 - T(i-1,i) - T(i+1,i);
    end
    T(2,1)     = kz(2)  * dt / (split*dz^2);
    T(1,1)     = 1.0 - T(2,1);
    T(nz-1,nz) = kz(nz) * dt / (split*dz^2);
    T(nz,nz)   = 1.0 - T(nz-1,nz);

    DIFF = ((T^split) * Cz - Cz) / dt;

end


return