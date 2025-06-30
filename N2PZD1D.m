more off
close all
clear all
format short g


% MATLAB / GNU OCTAVE -- PROGRAM "p20250510.m" 
%------------------------------------------------------------------------
% This script solves the basic N2PZD model in 1D (depth resolved 0 - 200m):
%
% dP/dt =  Fp - Ep - Gp - Mp; %Phytoplankton
% dZ/dt =  Fz - Ez - Gz - Mz; %Zooplankton
% dN/dt = -Fp + (epsPhy*Ep + epsZoo*Ez) + (omePhy*Mp + omeZoo*Mz) + md*D; %Nitrates
% dD/dt = (1-epsPhy)*Ep + (1-epsZoo)*Ez + (1-omePhy)*Mp + (1-omeZoo)*Mz - md*D; %Detritus
%
%------------------------------------------------------------------------
% AUTHORS 
% 
% Sergio M. Vallina (24 April 2022) - basic NPZD used in summer school
%
% Ziyu Guo and Katsumi Matsumoto, U. Minnesota (June, 2025) 
% (a) added Q10 temperature dependence to kinetic reaction rates
% (b) added vars: PO4, DIC, DOM (C, N, P) and D (P, C) w/o full C chemistry
% (c) linked phytoplankton C:N:P stoichiometry by the power law formula
% (d) added second P (phytoplankton); (possibly another Z) 
% (e) implemented Z grazing that depends on population and food quality
% (f) carry out experiments under different temperature regimes
%       - standing biomass, DOM and POM stoichiometry, nutrient cycling
%------------------------------------------------------------------------



%GLOBAL VARIABLES ========================================================
global knp kpp mup betap mp Isat epsPhy omePhy gzmax kgz betaz mz cz epsZoo omeZoo
global md mdom q10 npower H ZooCP_ref ZooCN_ref
global P2C0 sP2C_PO4 sP2C_NO3 sP2C_temp sP2C_I N2C0 sN2C_PO4 sN2C_NO3 sN2C_temp sN2C_I
global enable_pf sigma_pf tempK0 par0 PO40 NO30 kw kp wsink
global deltat C2P_uptake ntot0 ptot0 ctot0 stempC  
global zdepths ndepths deltaz mld parz0 KZ KZI PAR2D T jday jjday tresp
global INn IPn IZn IDn IDOMn INp IPp IZp IDp IDOMp INc IPc IZc IDc IDOMc IT
global Nn Pn Zn Dn DOMn Np Pp Zp Dp DOMp Nc Pc Zc Dc DOMc C2P_phyto
global Nndot Pndot Zndot Dndot DOMndot Npdot Ppdot Zpdot Dpdot DOMpdot Ncdot Pcdot Zcdot Dcdot DOMcdot
global pf_weighted1_t pf_weighted2_t C2P_uptake1_t C2P_uptake2_t C2P_phyto1_t C2P_phyto2_t C2P_detritus_t C2P_DOM_t C2P_zoo_t
global GPPn1_t GPPn2_t GSPn1_adj_t GSPn2_adj_t GSPn1_old_t GSPn2_old_t GSPn1_new_t GSPn2_new_t
global GPPc1_t GPPc2_t GSPc1_adj_t GSPc2_adj_t GSPc1_old_t GSPc2_old_t GSPc1_new_t GSPc2_new_t
global POCflux_t DOCflux_t
  

%MODEL KEYS ==============================================================

keyPhysics    = 'yes';  %If you want physical processes (diffusion and sinking).
% keyPhysics  = 'nop';  %If you do *not* want physical processes (diffusion and sinking).

keyModelDim   = '1D';   %Depth-resolved (several nodes in depth).
% keyModelDim = '0D';  %No depth-resolved (only one node).


%TEMPORAL RESOLUTION =====================================================
%Time step (dt) matters for the "ode4" solver because the user defines it.
%It does NOT matter for the "ode45" solver, which uses adaptative dt.

deltat  = 1/64;               %time step [days] for dz=5
%deltat  = 1/32;               %time step [days] for dz=10 
%deltat  = 1/8;               %time step [days] for dz=5 
t0      = deltat;
ndays   = 360;
years   = 5;
tmax    = ndays * years;
tspan   = [t0:deltat:tmax];
ntimes  = length(tspan);


%VERTICAL RESOLUTION =====================================================
deltaz = 5.0 %[m] (needs dt < 0.25)
%deltaz = 10.0 %[m] (needs dt < 0.25)

if strcmp(keyModelDim,'0D')
    zmin = deltaz;
    zmax = deltaz; %[m]
elseif strcmp(keyModelDim,'1D')
    zmin = 0.0;
    zmax = 200; %[m]
end

zdepths = [deltaz:deltaz:zmax]; %all nodes.
ndepths = length(zdepths); %number of vertical nodes.


%MODEL PARAMETERS ========================================================
% Photoplankton
knp    = [4.0  0.2];   % Nitrate half-saturation constant [mmolN*m-3]
kpp    = [0.4  0.04];  % Phosphate half-saturation constant [mmolP*m-3]
mup    = [3.0  0.5];   % Maximum uptake/growth rate [d-1]
betap  = [0.8  0.8];   % Assimilation efficiency [n.d.]
mp     = [0.05 0.05];  % Mortality specific rate [d-1]
Isat   = [120 120];    % Saturating irradiance [W*m-2]
epsPhy = [0.25 0.25];  % Exudation fraction going to nitrate [n.d.]
omePhy = [0.25 0.25];  % Mortality fraction going to nitrate [n.d.]

% Zooplanton
gzmax   = 1.0;     % Maximum ingestion rate [d-1]
kgz     = 0.75;    % Half-saturation ingestion constant [mmolN*m-3]
betaz   = 0.40;    % Assimilation efficiency [n.d.]
mz      = 0.05;    % Mortality rate [d-1]
cz      = 0.10;    % Closure mortality rate [d-1]
epsZoo  = 0.25;    % Exudation fraction going to nitrate [n.d.]
omeZoo  = 0.25;    % Mortality fraction going to nitrate [n.d.]

% Detritus / DOM / Enviormental Conditons
md     = 0.2;   % Detritus degradation rate [d-1]
mdom   = 0.2;   % DOM degradation rate [d⁻¹]
q10    = 2.0;   % Q10=2.0 (kinetic rates increase x2 for every 10C warming)
npower = 2.0;   % Grazing type: 1 = Holling type II; 2 = Holling Type III


% Redfield Ratios / Constants
RedCP         = 106 / 1;       % Redfield C:P
RedCN         = 106 / 16;      % Redfield C:N
RedNP         = 16 / 1;        % Redfield N:P
H             = 0.08;          % TT(2021)GRL 0.08; H=0 100% homeostatic; H=1 you're what you eat
ZooCP_ref     = RedCP;
ZooCN_ref     = RedCN;
umolkg2mmolm3 = (1e-3) / (1e-6) / 1025;  % unit conversion

% CNP power law constants (P1=eukaryotes, P2=cyanobacteria; from TM meta-analysis)
P2C0       = [11.6 6.3];       % permil: P/C=11.6 permil==>C:P=86.2; P/C=6.3 permil==>C:P=158.73
sP2C_PO4   = [0.58 0.28];
sP2C_NO3   = [0 0];
sP2C_temp  = [0 -8.0];
sP2C_I     = [0 0];

N2C0       = [151.0 151.0];    % permil: N/C=151.0 permil==>C:N=6.62
sN2C_PO4   = [0 0];
sN2C_NO3   = [0.22 0.22];
sN2C_temp  = [0 0];
sN2C_I     = [-0.05 -0.05];

% Grazing selectivity
enable_pf = 2;     % 0 = equal preference, fixed
                   % 1 = only prey density (squared) dependent (rho)
                   % 2 = only food quality-dependent (fq)
                   % 3 = both rho and fq dependent

sigma_pf  = 1.6;   % default = 1.0; option for when enable_pf=2 or 3
                   % selectiveness (smaller value peaks higher w/ narrower distribution)
                   % large sigma_pf --> enable_pf=1

% Environmental Conditions /Others
tempK0 = 291;                   % Reference temperature [K]
par0   = 70;                    % reference light level (W/m2)

PO40   = 0.57 * umolkg2mmolm3;  % reference PO4
NO30   = 5.70 * umolkg2mmolm3;  % reference NO3  

kw     = 0.04;                  % irradiance attenutation due to water [m-1]
kp     = 0.00;                  % irradiance attenutation due to phyto self-sheding [m2*mmolN-1]
%wsink  = 1.0;                   % [m*d-1]
wsink  = 1.0                   % [m*d-1]

% advection.m uses 1st order upwind differnecing...so check for assumption
% (e.g., dt=1/8 day and dz=10 m, wspeed could be up to 80 m/day)
if (wsink*deltat) > deltaz
    error('1st order upwind differencing assumption in particle sinking violated!')
end


%EXTERNAL FORCINGS =======================================================

%PHOTOSYNTHETIC ACTIVE RADIATION (PAR+) ----------------------------------
ttmax   = 360;                 % [days]
parmin  = 120;                 % winter PAR [W*m⁻²]
parmax  = 300;                 % summer PAR [W*m⁻²]

iparz0  = sincurve(parmin, parmax, ndays, ttmax); % Seasonal PAR
% iparz0 = parmax * ones(1, ttmax);               % Constant PAR version

%MIXED LAYER DEPTH (MLD) -------------------------------------------------
ttmax   = 360;                 %[days]
mldmin4sin = -30;                  %min MLD [m]
mldmin  = 30                  %min MLD [m]
mldmax  = 150;                 %max MLD [m]
[imld]  = sincurve(mldmax, mldmin4sin, ndays, ttmax);

%KM : create a summer plateau
mld_mask = imld > mldmin;
imld = imld.*mld_mask + mldmin*(1-mld_mask);

%SURFACE TURBULENCE DIFFUSION (SKZ)  -------------------------------------
kzmin      = 50;  %min surface turbulence [m^2*day-1]
kzmax      = 300; %max surface turbulence [m^2*day-1]
% kzmin      = 20;  %min surface turbulence [m^2*day-1]
% kzmax      = 200; %max surface turbulence [m^2*day-1]
kz_surface = kzmax - (kzmax - kzmin) * ((mldmax - imld) ./ (mldmax - mldmin)); %skz is a function of the MLD.

%TURBULENCE DIFFUSION PROFILES AT THE NODES j (KZ) -----------------------
kz_deep = ones(1, ndays);         %row-vect(1,ndays) daily minimum diffusivity [m2*day-1]
iKZ     = zeros(ndepths, ndays);

for iday = 1:ndays
    kz_profile    = diffprofile(kz_deep(iday), kz_surface(iday), imld(iday), zmin, zmax, deltaz);
                      
    iKZ(:, iday)  = kz_profile;   %Sigmoidal fitting to turbulent diffusivity profile -- [m2 * day-1]
end

%KZ AT MID-POINT BETWEEN NODES j-1/2 or j+1/2 (KZI) ----------------------
x       = 1:ndays;                                       %original time grid nodes.
z       = zdepths;                                       %original depth grid nodes.
xI      = 1:ndays;                                       %interpolated time grid nodes.
zI      = z(1) + 0.5*deltaz : deltaz : z(end) - 0.5*deltaz; %interpolated depth grid nodes.

[X , Z ]  = meshgrid(x , z );                            %convert to 2D original grids.
[XI, ZI]  = meshgrid(xI, zI);                            %convert to 2D interpolated grids.
iKZI      = interp2(X, Z, iKZ, XI, ZI);                  %linear interpolation of KZ values between nodes.

%MAKE LONGER ARRAYS (mld, par, KZ, KZI) ACCORDING TO NUMBER OF YEARS -----
KZ  = [];
KZI = [];
mld = [];
parz0 = []; 
for iyear = 1:years
    KZ  = [KZ, iKZ ];
    KZI = [KZI,iKZI];
    mld = [mld,imld];
    parz0 = [parz0,iparz0];
end
PAR2D = zeros(size(KZ));

%Surface Tempereture -----------------------------------------------------
ttmax = 360;
stempCmin = 12; %winter [deg C]; tresp=1 when temp 10C
stempCmax = 24;
[istempC] = sincurve(stempCmin,stempCmax,ndays,ttmax) ; %temperature (C) at the surface

for iyear = 1:years
    stempC = [stempC,istempC];
end


%STABILITY CONDITION FOR TURBULENCE EXPLICIT SCHEME ======================
%(IT DOES *NOT* APPLY WHEN USING AN IMPLICIT SCHEME)
KZmax     = max(KZ(:));
Stability = (deltaz^2) / (2*KZmax); % Condition for Numerical Stability of turbulent diffusion scheme. 
if (deltat*2) < Stability           % Check stability condition for turbulent diffusion.
    % disp('Okay: Condition of stability is alright')
else 
    disp('Error: Condition of stability violated!')
    disp('Solution: Please reduce dt or increase dz')
    return
end


%INTIAL CONDITIONS =======================================================
% Combine two types of phytoplankton
Pn0 = [0.75 0.75]; % Phytoplankton type 1,2 [mmolN*m-3]
Pp0 = Pn0 / RedNP;
Pc0 = Pp0 * RedNP;

% NPZD in Nitrogen units
Nn0   = 0.50; %Nutrient DIN [mmolN*m-3] KM: important for limitation and CNP
Zn0   = 0.25; %Zooplankton [mmolN*m-3]
Dn0   = 0.00; %Detrits PON [mmolN*m-3]
DOMn0 = 30/RedCN; %DON [mmolN*m-3];not Redfield, use Letscher's DOM stoichiometry to set

% NPZD in Phosphorus units - start out with Redfield relative to N
Np0   = Nn0/RedNP; %PO4 [mmolP*m-3]
Zp0   = Zn0/RedNP;
Dp0   = Dn0/RedNP; %Detritus P=POP
DOMp0 = 30/RedCP; %not Redfield, use Letscher's DOM stoichiometry to set

% NPZD in Carbon units
Nc0   = 2000; % DIC; mean surface ocean DIC~2000 umol/kg [mmolC/m3]
Zc0   = Zn0*RedCN;
Dc0   = Dn0*RedCN; %Detritus C=POC
DOMc0 = 30; % warm surface ocean labile/semilabile DOC umol/kg [mmolC/m3]

% Temperature
%T0 = 0;
T0 = stempCmin;

% ================
% Single-layer v0 
% ================
% Column vector of initial conditions for 1 depth layer (19)
v0 = [Nn0; Pn0'; Zn0; Dn0; DOMn0;
      Np0; Pp0'; Zp0; Dp0; DOMp0;
      Nc0; Pc0'; Zc0; Dc0; DOMc0;
      T0];
nvars = length(v0);  % Total number of state variables in 1-layer system

% ==============
%  Multilayer V0  
% ==============
vdepths = ones(ndepths, 1);  % depth profile multiplier

% for nutrient n/p/c: (ndepths × 6; matrix multiplication; 20x6 each=360)
dn = vdepths * [Nn0, Pn0, Zn0, Dn0, DOMn0];
dp = vdepths * [Np0, Pp0, Zp0, Dp0, DOMp0];
dc = vdepths * [Nc0, Pc0, Zc0, Dc0, DOMc0];

% for temperature (ndepths x 1; 20x1)
dT0 = vdepths * T0;

% Flatten and concatenate all to form the column vector V0 (360)
V0 = [dn(:); dp(:); dc(:); dT0];

ntot0 = sum(V0(1:6*ndepths));   %initial total N mass [mmolN*m-3]   % total mass for n,p,c
ptot0 = sum(V0(6*ndepths+1:12*ndepths));  %initial total P mass [mmolN*m-3]
ctot0 = sum(V0(12*ndepths+1:18*ndepths)); %initial total C mass [mmolN*m-3]

% Variable indices
INn = 1:ndepths;                    % nutrient N
IPn = (ndepths)+1:(ndepths*3); 
IZn = (ndepths*3)+1:(ndepths*4);    
IDn = (ndepths*4)+1:(ndepths*5);
IDOMn = (ndepths*5)+1:(ndepths*6);

INp = (ndepths*6)+1:(ndepths*7);    % nutrient P
IPp = (ndepths*7)+1:(ndepths*9); 
IZp = (ndepths*9)+1:(ndepths*10);
IDp = (ndepths*10)+1:(ndepths*11); 
IDOMp = (ndepths*11)+1:(ndepths*12);

INc = (ndepths*12)+1:(ndepths*13);  % nutrient C
IPc = (ndepths*13)+1:(ndepths*15);
IZc = (ndepths*15)+1:(ndepths*16); 
IDc = (ndepths*16)+1:(ndepths*17);
IDOMc = (ndepths*17)+1:(ndepths*18); 

IT = (ndepths*18)+1:(ndepths*19);   % temperature


%IF YOU WANT TO REMOVE ALL PHYSICAL PROCESSES ----------------------------
if strcmp(keyPhysics,'nop')
    KZ (:,:) = 1;
    KZI(:,:) = 1;
    wsink = 0;
end


%INITIALIZE TRANSIENT OUPUT VARIABLES...every day ------------------------
pf_weighted1_t  = zeros(ndepths,tmax);
pf_weighted2_t  = zeros(ndepths,tmax);

C2P_uptake1_t   = zeros(ndepths,tmax);
C2P_uptake2_t   = zeros(ndepths,tmax);
C2P_phyto1_t    = zeros(ndepths,tmax);
C2P_phyto2_t    = zeros(ndepths,tmax);
C2P_zoo_t       = zeros(ndepths,tmax);
C2P_detritus_t  = zeros(ndepths,tmax);
C2P_DOM_t       = zeros(ndepths,tmax);

GPPn1_t         = zeros(ndepths,tmax);
GPPn2_t         = zeros(ndepths,tmax);
GSPn1_old_t     = zeros(ndepths,tmax);
GSPn2_old_t     = zeros(ndepths,tmax);
GSPn1_new_t     = zeros(ndepths,tmax);
GSPn2_new_t     = zeros(ndepths,tmax);
GSPn1_adj_t     = zeros(ndepths,tmax);
GSPn2_adj_t     = zeros(ndepths,tmax);

GPPc1_t         = zeros(ndepths,tmax);
GPPc2_t         = zeros(ndepths,tmax);
GSPc1_old_t     = zeros(ndepths,tmax);
GSPc2_old_t     = zeros(ndepths,tmax);
GSPc1_new_t     = zeros(ndepths,tmax);
GSPc2_new_t     = zeros(ndepths,tmax);
GSPc1_adj_t     = zeros(ndepths,tmax);
GSPc2_adj_t     = zeros(ndepths,tmax);

POCflux_t       = zeros(ndepths,tmax);
DOCflux_t       = zeros(ndepths,tmax);


%Experiments =============================================================
% enable_pf_values = [0];   % single runs
% sigma_pf_values = [0];

enable_pf_values = [0   1   2   2   2   2  ]; 
sigma_pf_values  = [0.0 0.0 0.7 1.0 1.3 2.0];

numExperiments = numel(enable_pf_values); 

plotData = struct();    % create a struct to store output from all the experiments

for i = 1 : numExperiments

    enable_pf = enable_pf_values(i);
    sigma_pf =  sigma_pf_values(i);
    ['experiment i=',num2str(i),', pf=',num2str(enable_pf),', sigma=' num2str(sigma_pf)]

    %SOLVE EQUATIONS USING RUNGE-KUTTA ODE45 =============================
    %---------------------------------------------------------------------
    % NormControl --> Requests ode45 to control the error in each step, so that 
    % the norm of the error is smaller than Reltol*norm(solution) or AbsTol, 
    % whichever is largest; norm(e) <= max( RelTol*norm(y), AbsTol ). 
    % By default ode45 uses a more stringent componentwise error control.
    %---------------------------------------------------------------------
    ode45options = odeset('AbsTol',1e-12,'RelTol',1e-6,'NonNegative',[1:ndepths*nvars]); %Slower but better
    %%ode45options = odeset ('AbsTol',1e-12,'RelTol',1e-6,'NonNegative',[1:ndepths*nvars],'NormControl','on'); %Faster but worse
    tic 
    jday  = 1;
    jjday = 1;
 
    [Vode] = ode4(@N2PZD1D_eqs,tspan,V0);
    % [tode,Vode] = ode45(@NPZD1D_eqs,tspan,V0,ode45options);
    Vode = Vode'; %transpose matrix [rows-StateVariables, columns-time]
    Vode_reshape = reshape(Vode,[ndepths,nvars,ntimes]);
    
    %GET DAILY AVERAGES ==================================================
    [aveVode,avetspan] = dailyAverages(Vode,deltat);
    deltaday = 1; %[day] 

    %STATE VARIABLES =====================================================
    tode = avetspan;   % time                  

    Nnode = aveVode(INn,:);     % for nutrient N
    Pnode = aveVode(IPn,:);     % 2 types of phytoplankton
    Znode = aveVode(IZn,:); 
    Dnode = aveVode(IDn,:);  
    DOMnode = aveVode(IDOMn,:); 

    Npode = aveVode(INp,:);     % for nutrient P
    Ppode = aveVode(IPp,:); 
    Zpode = aveVode(IZp,:); 
    Dpode = aveVode(IDp,:); 
    DOMpode = aveVode(IDOMp,:); 

    Ncode = aveVode(INc,:);     % for nutrient C
    Pcode = aveVode(IPc,:); 
    Zcode = aveVode(IZc,:); 
    Dcode = aveVode(IDc,:); 
    DOMcode = aveVode(IDOMc,:); 

    Tode = aveVode(IT,:);       % temperature (set up the boundary)
    
    J = [(ndays/deltaday)*(years-1)+1:(ndays/deltaday)*years]; %Last year days (ie, days 1441~1800)
    
    % save output struct
    plotData(i).NnodeSSP = Nnode(:,J); 
    plotData(i).PnodeSSP = Pnode(:,J);      % twice as long; 2 P types
    plotData(i).ZnodeSSP = Znode(:,J); 
    plotData(i).DnodeSSP = Dnode(:,J); 
    plotData(i).DOMnodeSSP = DOMnode(:,J);
    
    plotData(i).NpodeSSP = Npode(:,J); 
    plotData(i).PpodeSSP = Ppode(:,J); 
    plotData(i).ZpodeSSP = Zpode(:,J); 
    plotData(i).DpodeSSP = Dpode(:,J); 
    plotData(i).DOMpodeSSP = DOMpode(:,J); 
    
    plotData(i).NcodeSSP = Ncode(:,J);  
    plotData(i).PcodeSSP = Pcode(:,J); 
    plotData(i).ZcodeSSP = Zcode(:,J); 
    plotData(i).DcodeSSP = Dcode(:,J); 
    plotData(i).DOMcodeSSP = DOMcode(:,J); 
    
    plotData(i).TodeSSP = Tode(:,J); 
    plotData(i).mldSSP = mld(:,J);
    plotData(i).kzSSP = iKZ;
    plotData(i).parSSP = PAR2D(:,J);
    plotData(i).sigma = sigma_pf;

    plotData(i).pf_weighted1_tSSP = pf_weighted1_t(:,J); 
    plotData(i).pf_weighted2_tSSP = pf_weighted2_t(:,J);
    plotData(i).C2P_phyto1_tSSP = C2P_phyto1_t(:,J); 
    plotData(i).C2P_phyto2_tSSP = C2P_phyto2_t(:,J);
    plotData(i).C2P_zoo_tSSP = C2P_zoo_t(:,J); 
    plotData(i).C2P_detritus_tSSP = C2P_detritus_t(:,J); 
    plotData(i).C2P_DOM_tSSP = C2P_DOM_t(:,J);

    plotData(i).GPPn1_tSSP = GPPn1_t(:,J);
    plotData(i).GPPn2_tSSP = GPPn2_t(:,J);
    plotData(i).GSPn1_old_tSSP = GSPn1_old_t(:,J); 
    plotData(i).GSPn2_old_tSSP = GSPn2_old_t(:,J); 
    plotData(i).GSPn1_new_tSSP = GSPn1_new_t(:,J); 
    plotData(i).GSPn2_new_tSSP = GSPn2_new_t(:,J); 
    plotData(i).GSPn1_adj_tSSP = GSPn1_adj_t(:,J);
    plotData(i).GSPn2_adj_tSSP = GSPn2_adj_t(:,J);

    plotData(i).GPPc1_tSSP = GPPc1_t(:,J);
    plotData(i).GPPc2_tSSP = GPPc2_t(:,J);
    plotData(i).GSPc1_old_tSSP = GSPc1_old_t(:,J); 
    plotData(i).GSPc2_old_tSSP = GSPc2_old_t(:,J); 
    plotData(i).GSPc1_new_tSSP = GSPc1_new_t(:,J); 
    plotData(i).GSPc2_new_tSSP = GSPc2_new_t(:,J); 
    plotData(i).GSPc1_adj_tSSP = GSPc1_adj_t(:,J);
    plotData(i).GSPc2_adj_tSSP = GSPc2_adj_t(:,J);

    plotData(i).POCfluxSSP = POCflux_t(:,J); 
    plotData(i).DOCfluxSSP = DOCflux_t(:,J);

end

% SAVE and PLOT ==========================================================


% filename=['1d_dz' num2str(deltaz) '_mld' num2str(mldmin) '_w' num2str(wsink) '.mat']
% eval(['save ',filename...
%     ' plotData ndays deltaday zmax deltaz numExperiments ndepths'])

return

%run('plot_all.m')
%run('plot_all_km.m')
run('plot_km.m')
