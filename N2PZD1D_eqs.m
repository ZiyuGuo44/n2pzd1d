function [Vdot]=N2PZD1D_eqs(iTime,V0)

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
  

%STATE VARIABLES =========================================================
Nn    = V0(INn);  %NO3 [mmolN*m-3]
Pn    = V0(IPn);
Zn    = V0(IZn);  %Zooplankton [mmolN*m-3]
Dn    = V0(IDn);  %Detritus [mmolN*m-3]
DOMn  = V0(IDOMn);  %DON

Np    = V0(INp);  %PO4 [mmolP*m-3]
Pp    = V0(IPp);
Zp    = V0(IZp); %Z
Dp    = V0(IDp); %D
DOMp  = V0(IDOMp); %DOP

Nc    = V0(INc); %DIC [mmolC*m-3]
Pc    = V0(IPc);
Zc    = V0(IZc); %Zooplankton Carbon content
Dc    = V0(IDc); %Detritus Carbon content
DOMc  = V0(IDOMc); %DOC
T     = V0(IT);

Pn = reshape(V0(IPn), ndepths, 2); 
Pp = reshape(V0(IPp), ndepths, 2); 
Pc = reshape(V0(IPc), ndepths, 2); 


%DAY OF SIMULATION =======================================================
jday0 = jday; %Previous day.
jday  = ceil(iTime); %Current day.
jday1 = jday; %Current day.

if jday0 == jday1
    newday = 'nop';
    jjday  = jday0; %Previous day.
else
    newday = 'yes';
    jjday  = jday1; %Curent day.
end


%CHECK MASS CONSERVATION & NEGATIVE VALUES ===============================
ntoti = sum(V0(1:6*ndepths));  
ptoti = sum(V0(6*ndepths+1:12*ndepths));  
ctoti = sum(V0(12*ndepths+1:18*ndepths)); 

ndist = ntoti - ntot0;
pdist = ptoti - ptot0;
cdist = ctoti - ctot0;
ndistmax = 1e-6;

if abs(ndist) > ndistmax
    masscheck = [iTime,ntot0,ntoti,ndist];
    disp(masscheck);
    disp('Error!!! N mass is NOT conserved!')
    if abs(pdist) > ndistmax
        masscheck = [iTime,ptot0,ptoti,ndist];
        disp('Error!!! P mass is NOT conserved!')
        if abs(cdist) > ndistmax
            masscheck = [iTime,ctot0,ctoti,ndist];
            disp('Error!!! C mass is NOT conserved!')
            pause
        end
        pause
    end
    pause
end

Ineg = find(V0(:) < 0);
if length(Ineg) > 0
    disp(iTime);
    disp('Error!!! there are NEGATIVE concentrations!')
    pause
end


%TURBULENT DIFFUSION =====================================================
kz  = KZ (:,jday);
kzI = KZI(:,jday);

% keyScheme = 'Explicit';
keyScheme = 'Implicit';

DIFF_Nn   = diffusion_i(Nn,   deltat, deltaz, kz, ndepths, keyScheme);
DIFF_Pn   = diffusion_i(Pn,   deltat, deltaz, kz, ndepths, keyScheme);
DIFF_Zn   = diffusion_i(Zn,   deltat, deltaz, kz, ndepths, keyScheme);
DIFF_Dn   = diffusion_i(Dn,   deltat, deltaz, kz, ndepths, keyScheme);
DIFF_DOMn = diffusion_i(DOMn, deltat, deltaz, kz, ndepths, keyScheme);

DIFF_Np   = diffusion_i(Np,   deltat, deltaz, kz, ndepths, keyScheme);
DIFF_Pp   = diffusion_i(Pp,   deltat, deltaz, kz, ndepths, keyScheme);
DIFF_Zp   = diffusion_i(Zp,   deltat, deltaz, kz, ndepths, keyScheme);
DIFF_Dp   = diffusion_i(Dp,   deltat, deltaz, kz, ndepths, keyScheme);
DIFF_DOMp = diffusion_i(DOMp, deltat, deltaz, kz, ndepths, keyScheme);

DIFF_Nc   = diffusion_i(Nc,   deltat, deltaz, kz, ndepths, keyScheme);
DIFF_Pc   = diffusion_i(Pc,   deltat, deltaz, kz, ndepths, keyScheme);
DIFF_Zc   = diffusion_i(Zc,   deltat, deltaz, kz, ndepths, keyScheme);
DIFF_Dc   = diffusion_i(Dc,   deltat, deltaz, kz, ndepths, keyScheme);
DIFF_DOMc = diffusion_i(DOMc, deltat, deltaz, kz, ndepths, keyScheme);

DIFF_T    = diffusion_i(T,    deltat, deltaz, kz, ndepths, keyScheme);


%VERTICAL SINKING ========================================================
ADV_Nn   = 0;
ADV_Pn   = 0;
ADV_Zn   = 0;
ADV_Dn   = advection(Dn, deltat, deltaz, wsink, ndepths);
ADV_DOMn = 0;

ADV_Np   = 0;
ADV_Pp   = 0;
ADV_Zp   = 0;
ADV_Dp   = advection(Dp, deltat, deltaz, wsink, ndepths);
ADV_DOMp = 0;

ADV_Nc   = 0;
ADV_Pc   = 0;
ADV_Zc   = 0;
ADV_Dc   = advection(Dc, deltat, deltaz, wsink, ndepths);
ADV_DOMc = 0;

ADV_T   = 0;


%TEMPERATURE =============================================================
sst  = stempC(1,jday);
dsst = sst - T(1);
T(1) = sst;

Tdot = zeros(ndepths,1);
Tdot =  Tdot + DIFF_T + ADV_T;

Tdot(1) = dsst;


%MODEL TERMS =============================================================
jpar0  = parz0(jday); %Photo. Active. Radiation at the surface [W*m-2]
jPAR   = jpar0*exp(-kw*zdepths(:)); %[W*m-2] PAR profile.
jTEMPC = T ;  % Temp 
jtempK = jTEMPC + 273.15; % Kelvin
tresp  = q10.^((jTEMPC(:)-10)/10)-q10.^((jTEMPC(:)-32)/3); %Temp response after Blackford 2004


%Phytoplankton Equations =================================================
% Phytoplankton limitation 
Qpar = jPAR./ Isat.*exp(1 - jPAR./Isat); %Phy light limitation [n.d.] values between 0 and 1.
Qdin = Nn ./ (knp + Nn);%Phytoplankton NO3 limitation [n.d.] values between 0 and 1
Qdip = Np ./ (kpp + Np); %Phytoplankton PO4 limitation [n.d.] values between 0 and 1

% Phytoplankton uptake C:N:P ratio according to TT's meta-analysis; permil
P2Cpermil = P2C0 .* (Np/PO40).^sP2C_PO4.*(Nn/NO30).^sP2C_NO3.*(jtempK/tempK0).^sP2C_temp.*(jPAR./par0).^sP2C_I;
N2Cpermil  = N2C0  .* (Np/PO40).^sN2C_PO4.*(Nn/NO30).^sN2C_NO3.*(jtempK/tempK0).^sN2C_temp.*(jPAR/par0).^sN2C_I;
C2P = 1./(P2Cpermil/1000);
C2N = 1./(N2Cpermil/1000);

C2P_uptake = min(546.7,max(C2P,26.6));
C2N_uptake = min(30,max(C2N,2));
N2P_uptake = C2P_uptake ./ C2N_uptake;

%Phytoplankton cellular quota C:N:P ratio
C2P_phyto = Pc./Pp;
C2N_phyto = Pc./Pn;
N2P_phyto = Pn./Pp;
C2P_uptake1_t(:,jday) = C2P_uptake(:,1);
C2P_uptake2_t(:,jday) = C2P_uptake(:,2);
C2P_phyto1_t(:,jday) = C2P_phyto(:,1);
C2P_phyto2_t(:,jday) = C2P_phyto(:,2);

is_N_limit = (Nn .* Qdin) <= (Np .* Qdip .* N2P_uptake);
Qnut = Qdin;
Qnut(~is_N_limit) = Qdip(~is_N_limit);

GPPn = (mup.* tresp).*Qpar.*Qnut.*Pn;    %Phy gross primary production
EPn  = (1- betap).*GPPn;   %Phy exudation that increases w/ temperature...5-20%; betap is assimilated
MPn  = (mp.* tresp).*Pn; 
GPPp = GPPn./N2P_uptake;                 %Linked to P by N:P uptake stoichiometry
EPp  = (1- betap).*GPPp; 
MPp  = MPn./N2P_phyto;                   %Linked to P by phytoplankton N:P              
GPPc = GPPn.*C2N_uptake;                 %Linked to C by C:N uptake stoichiometry
EPc  = (1- betap).*GPPc;
MPc  = MPn.*C2N_phyto;                   %Linked to P by phytoplankton C:P                            


%Zooplankton Equations ===================================================
C2N_zoo = Zc./Zn;                                         % Zoo cellular quota C:N:P ratio
N2P_zoo = Zn./Zp;
C2P_zoo = C2N_zoo .* N2P_zoo;
Qphyn   = Pn.^npower ./ (kgz^npower + sum(Pn.^npower,2)); % Zoo grazing limitation [n.d.] for P>1
C2P_zoo_t(:,jday) = C2P_zoo;

% (1) Various grazing preferences ----------------------------------------
%1) Nominal preference (equal=no preference)            
pf_nominal = 0.5*ones(size(Qphyn)); 

%2a) Density (linear) dependent grazing (Fasham; denominator is total food)
% e.g., pf_nominal(1)*Pn(1) / (pf_nominal(1)*Pn(1) + pf_nominal(2)*Pn(2))
pf_rho1 = pf_nominal.*Pn./(sum(pf_nominal.*Pn,2));

%2b) Density (squared) dependent grazing preference (as in NPPZ0D model)
% e.g., Pn(1)^2 ./  (Pn(1)^2 + Pn(2)^2); assumes pf_nominal=[0.5 0.5]
% density-squared dependence comes from independent encounter and interaction probabilities
pf_rho2 = (pf_nominal.*Pn).^2 ./  sum((pf_nominal.*Pn).^2,2);

%3) Food quality dependent grazing preference
C2P_optimal = C2N_zoo .* N2P_zoo; % optimal CP is zooplankton C:P, assuming gross growth efficiency 
%C2P_optimal = 106;
x  = (C2P_phyto - C2P_optimal) ./ C2P_optimal;
fq = normpdf(x,0.0,sigma_pf);  

if enable_pf == 0
    pf_weighted = pf_nominal; % sum to 1
elseif enable_pf == 1
    pf_weighted = pf_rho2;
elseif enable_pf == 2
    pf_weighted = fq;
elseif enable_pf == 3
    pf_weighted = fq .* pf_rho2;
end
pf_weighted = pf_weighted ./ sum(pf_weighted,2);
pf_weighted1_t(:,jday) = pf_weighted(:,1);
pf_weighted2_t(:,jday) = pf_weighted(:,2);

% (2) Equations ----------------------------------------------------------
GSPn  = (gzmax.*tresp).*pf_weighted.*Qphyn.* Zn; %Zoo grazing (gross secondary production) [mmolN*m-3*d-1]
MZn   = (mz.*tresp).*Zn;                         %Zoo natural mortality [mmolN*m-3*d-1]
GTPn  = (cz.*tresp).*(Zn.*Zn);                   %Higher order predatation on Zoo (gross tertiary prod; closure) [mmolN*m-3*d-1]

GSPp  = GSPn./N2P_phyto;                         %GSPp is related to GSPn by phytolankton N:P ratio
MZp   = MZn./N2P_zoo;                            %Linked to P by zooplankton N:P ratio
GTPp  = GTPn./N2P_zoo;                           %Linked to P by zooplankton N:P ratio

GSPc  = GSPn.*C2N_phyto;                         %GSPc is related to GSPn by cellular C:N ratio
MZc   = MZn.*C2N_zoo;                            %Linked to P by zooplankton N:P ratio
GTPc  = GTPn.*C2N_zoo;  

% (3) Zooplankton homeostatic adjustment of GSP fluxes -------------------
% First, adjust either N or P whichever is more unbalanced from ideal zoo C:N and C:P
P2C_zoo_ideal = (1/ZooCP_ref)^(1-H) .* (1./C2P_phyto).^H;
N2C_zoo_ideal = (1/ZooCN_ref)^(1-H) .* (1./C2N_phyto).^H;

Nrelease_byZ = max(0, GSPn - GSPc.*N2C_zoo_ideal);
Prelease_byZ = max(0, GSPp - GSPc.*P2C_zoo_ideal);

Nlimit_vsP_1 = GSPn - (N2C_zoo_ideal./P2C_zoo_ideal).*GSPp; %N excess if >0
Plimit_vsN_1 = GSPp - (P2C_zoo_ideal./N2C_zoo_ideal).*GSPn; %P excess if >0

Nadj_flux_1 = Nrelease_byZ.*double(Nlimit_vsP_1>0); %Either N release
Padj_flux_1 = Prelease_byZ.*double(Plimit_vsN_1>0); % or P release, not both

GSPn_1 = GSPn - Nadj_flux_1; %Net change after first adjustment
GSPp_1 = GSPp - Padj_flux_1;

% Second, adjust C flux relative to N or P whichever is more limiting
Cflux_vsNP  = max(GSPc-GSPn_1./N2C_zoo_ideal, GSPc-GSPp_1./P2C_zoo_ideal);
Cadj_flux_1 = Cflux_vsNP.*double(Cflux_vsNP>0);

% Third, go back to N and P, and adjust N:P flux
Nlimit_vsP_2 = GSPn_1 - (N2C_zoo_ideal./P2C_zoo_ideal).*GSPp_1; %N excess: P limitation if >0
Plimit_vsN_2 = GSPp_1 - (P2C_zoo_ideal./N2C_zoo_ideal).*GSPn_1; %P excess: N limitation if >0

Nadj_flux_2 = Nlimit_vsP_2.*double(Nlimit_vsP_2>0);
Padj_flux_2 = Plimit_vsN_2.*double(Plimit_vsN_2>0);

GSPn_adj = Nadj_flux_1 + Nadj_flux_2; %Total homeostatic adjustment in flux (not biomass)
GSPp_adj = Padj_flux_1 + Padj_flux_2;
GSPc_adj = Cadj_flux_1;

GSPn_2 = nanmax(0,GSPn - GSPn_adj); %New GSP after homeostatic adjustment; 
GSPp_2 = nanmax(0,GSPp - GSPp_adj);
GSPc_2 = nanmax(0,GSPc - GSPc_adj);

% KM change output
GPPc1_t(:,jday)     = GPPc(:,1);        % GPP of species 1 based on C 
GPPc2_t(:,jday)     = GPPc(:,2); 
GSPc1_adj_t(:,jday) = GSPc_adj(:,1);    % Homeostatic adj flux
GSPc2_adj_t(:,jday) = GSPc_adj(:,2);
GSPc1_old_t(:,jday) = GSPc(:,1);        % Old GSP before homeo adj
GSPc2_old_t(:,jday) = GSPc(:,2);
GSPc1_new_t(:,jday) = GSPc_2(:,1);      % New GSP after homeo adj
GSPc2_new_t(:,jday) = GSPc_2(:,2);

GPPn1_t(:,jday)     = GPPn(:,1);        % GPP of species 1 based on N 
GPPn2_t(:,jday)     = GPPn(:,2); 
GSPn1_adj_t(:,jday) = GSPn_adj(:,1);
GSPn2_adj_t(:,jday) = GSPn_adj(:,2);
GSPn1_old_t(:,jday) = GSPn(:,1);
GSPn2_old_t(:,jday) = GSPn(:,2);
GSPn1_new_t(:,jday) = GSPn_2(:,1);
GSPn2_new_t(:,jday) = GSPn_2(:,2);

EZn = (1- betaz).*GSPn_2; %Z exudation calculated on the basis of adjusted GSP
EZp = (1- betaz).*GSPp_2;
EZc = (1- betaz).*GSPc_2;


%Detritus Equations ======================================================
DDn = (md*tresp).*Dn;                 %PON degradation rate;
DDp = (md*tresp).*Dp;                 %POP degradation rate;
DDc = (md*tresp).*Dc;                 %POC degradation rate;
C2P_detritus = Dc./Dp;
C2P_detritus_t(:,jday) = C2P_detritus;


%DOM Equations ===========================================================
DDOMn = (mdom*tresp).*DOMn;          %DON degradation rate;
DDOMp = (mdom*tresp).*DOMp;          %DOP degradation rate;
DDOMc = (mdom*tresp).*DOMc;          %DOC degradation rate;
C2P_DOM = DOMc./DOMp;
C2P_DOM_t(:,jday) = C2P_DOM;


%ECOSYSTEM MODEL EQUATIONS ===============================================
% Nitrogen NPZD + DOM [mmolN*m-3*d-1]
Nndot   = DDOMn - sum(GPPn,2);
Pndot   = GPPn - EPn - GSPn - MPn; 
Zndot   = sum(GSPn_2,2) - sum(EZn,2) - GTPn - MZn; 
Dndot   = GTPn + sum((1-epsPhy).* EPn,2) + (1-epsZoo).* sum(EZn + GSPn - GSPn_2,2) + sum((1-omePhy).*MPn,2) + (1-omeZoo).*MZn - DDn; 
DOMndot = sum(epsPhy.* EPn,2) + epsZoo.* sum(EZn + GSPn - GSPn_2,2) +  sum(omePhy.* MPn,2) + omeZoo.*MZn + DDn - DDOMn;

% Phosphorus NPZD + DOM [mmolP*m-3*d-1]
Npdot   = DDOMp - sum(GPPp,2); 
Ppdot   = GPPp - EPp - GSPp - MPp; 
Zpdot   = sum(GSPp_2,2) - sum(EZp,2) - GTPp - MZp; 
Dpdot   = GTPp + sum((1-epsPhy).* EPp,2) + (1-epsZoo).* sum(EZp + GSPp - GSPp_2,2) + sum((1-omePhy).*MPp,2) + (1-omeZoo).*MZp - DDp; 
DOMpdot = sum(epsPhy.* EPp,2) + epsZoo.* sum(EZp + GSPp - GSPp_2,2) +  sum(omePhy.* MPp,2) + omeZoo.*MZp + DDp - DDOMp;

% Carbon NPZD + DOM [mmolC*m-3*d-1]
Ncdot   = DDOMc - sum(GPPc,2); 
Pcdot   = GPPc - EPc - GSPc - MPc; 
Zcdot   = sum(GSPc_2,2) - sum(EZc,2) - GTPc - MZc; 
Dcdot   = GTPc + sum((1-epsPhy).* EPc,2) + (1-epsZoo).* sum(EZc + GSPc - GSPc_2,2) + sum((1-omePhy).*MPc,2) + (1-omeZoo).*MZc - DDc; 
DOMcdot = sum(epsPhy.* EPc,2) + epsZoo.* sum(EZc + GSPc - GSPc_2,2) +  sum(omePhy.* MPc,2) + omeZoo.*MZc + DDc - DDOMc;


%CHECK MASS CONSERVATION =================================================
sumODEs  = Nndot+  sum(Pndot,2) + Zndot + Dndot + DOMndot +  Npdot +  sum(Ppdot,2) + Zpdot + Dpdot + DOMpdot + ...
           Ncdot +  sum(Pcdot,2) + Zcdot + Dcdot + DOMcdot;
xdist    = abs(sumODEs - 0d0);
xdistmax = 1d-4;
I = find(xdist > xdistmax);

if length(I) > 1
    sumODEs %the sum of all ODEs should be equal to zero (closed system)!
    disp('Error!!!! Mass is NOT conserved!!!')
    pause
end


%ADD DIFFUSION + ADVECTION ===============================================
Nndot    = Nndot    + DIFF_Nn   + ADV_Nn;
Pndot    = Pndot    + DIFF_Pn   + ADV_Pn;
Zndot    = Zndot    + DIFF_Zn   + ADV_Zn;
Dndot    = Dndot    + DIFF_Dn   + ADV_Dn;
DOMndot  = DOMndot  + DIFF_DOMn + ADV_DOMn;

Npdot    = Npdot    + DIFF_Np   + ADV_Np;
Ppdot    = Ppdot    + DIFF_Pp   + ADV_Pp;
Zpdot    = Zpdot    + DIFF_Zp   + ADV_Zp;
Dpdot    = Dpdot    + DIFF_Dp   + ADV_Dp;
DOMpdot  = DOMpdot  + DIFF_DOMp + ADV_DOMp;

Ncdot    = Ncdot    + DIFF_Nc   + ADV_Nc;
Pcdot    = Pcdot    + DIFF_Pc   + ADV_Pc;
Zcdot    = Zcdot    + DIFF_Zc   + ADV_Zc;
Dcdot    = Dcdot    + DIFF_Dc   + ADV_Dc;
DOMcdot  = DOMcdot  + DIFF_DOMc + ADV_DOMc;


% OrgC FLUX ==============================================================
% advection.m is based on first order upwind method, which assumes that 
% particles only enter from the immediate upstream cell: wspeed x dt < dz

POCflux = Dc * wsink;   % mmolC/m2/d (instantaneous flux OUT OF this box)

D           = kz;                       % m^2*day-1
dc          = zeros(size(DOMc));
dc(2:end-1) = (DOMc(3:end) - DOMc(1:end-2)) / 2;
dc(1)       =  DOMc(2)     - DOMc(1);
dc(end)     =  DOMc(end)   - DOMc(end-1);

DOCflux = -D .* dc ./ deltaz;       % mmolC/m2/d (instantaneous flux)          

POCflux_t(:, jday) = POCflux;       % mmolC/m2/day
DOCflux_t(:, jday) = DOCflux;       % mmolC/m2/day


%OUTPUT ==================================================================
PAR2D(:,jday) = jPAR; %output variable for ploting.

Vdot = [Nndot;Pndot(:);Zndot;Dndot;DOMndot; 
        Npdot;Ppdot(:);Zpdot;Dpdot;DOMpdot;
        Ncdot;Pcdot(:);Zcdot;Dcdot;DOMcdot;Tdot]; % Must be a *column* vector


return