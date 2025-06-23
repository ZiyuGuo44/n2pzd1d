function [DIFF] = diffusion_explicit(C0,dt,dz,kz,kzI,zjmax)

C0  = C0(:)'; %row vector (0-Zm).
kz  = kz(:)'; %row vector (0-Zm).
kzI = kzI(:)'; %row vector (0-Zm).

ScaleFactor = (1 /dz^2); %Okay 
%%ScaleFactor = (dt/dz^2); %Wrong!!!! (do NOT use) 

%%%%%%%%%%%%%%%%%%%%%%%
%FOR THE SURFACE (j=1):
%%%%%%%%%%%%%%%%%%%%%%%
J = 1;
%C0surface = 0; %Frontier Conition (absorbent frontier)
C0surface = C0(J); %Frontier Conition (reflectance frontier)
kzsurface = kz(J);
DIFFsurface = (kzI(J)*(C0(J+1)-C0(J)) - kzsurface*(C0(J)-C0surface)) * ScaleFactor; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MIDDLE DEPTHS (j=2:zjmax-1):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = [2:zjmax-1];
DIFFmiddle = (kzI(J).*(C0(J+1)-C0(J)) - kzI(J-1).*(C0(J)-C0(J-1))) * ScaleFactor;

%%%%%%%%%%%%%%%%%%
%BOTTOM (j=zjmax):
%%%%%%%%%%%%%%%%%%
J = [zjmax];
%C0deep = 0; %C.F. (absorb)
C0deep = C0(J); %F.C. (reflect)
kzdeep = kz(J);
DIFFdeep = (kzdeep*(C0deep-C0(J)) - kzI(J-1)*(C0(J)-C0(J-1))) * ScaleFactor;

%%%%%%%%%%%%%%%
%TOTAL PROFILE:
%%%%%%%%%%%%%%%
DIFF = [DIFFsurface,DIFFmiddle,DIFFdeep];
DIFF = DIFF(:); %column vector.


return