%FHFUNCTIONS
% This routines generates a set of functions to compute flash heating
% parameters. Require the following files:
%    - Thermal_properties_Qz_Olv.txt (from Clauser and Huenges, 1995)
%    - Thermal_properties_rocks.txt (from Vosteen and SchellSchmidt, 2003)
%    - Hardness_Olivine_T.txt (from Evans and Goetze, 1979)
%    - Hardness_Qz_T.txt (from Evans, 1984)
%
%and makes use of the function flashheatingsoap_G.
%
% Output:
%   kT_qz:   thermal conductivity of qtz, function of temperature.
%   ath_qz:  thermal diffusivity of qtz, function of temperature.
%   rhoc_qz: heat capacity of qtz, function of temperature.
%   kT_olv:  thermal conductivity of olivine, function of temperature.
%   ath_olv: thermal diffusivity of olivine, function of temperature.
%   rhoc_olv:  heat capacity of olivine, function of temperature.
%   kT_rock:   thermal conductivity of bulk rock, function of temperature.
%   ath_rock:  thermal diffusivity of bulk rock, function of temperature.
%   rhoc_rock: heat capacity of bulk rock, function of temperature.
%   sy_qz:   yield strength of quartz, function of temperature.
%   sy_olv:  yield strength of olivine, function of temperature.
%   Vw:      weakening velocity. Function of heat diffusivity (ath),
%            asperity size (D), heat capacity (rhoc), weakening temperature
%            (Tw), initial temperature (T0), asperity shear strength
%            (tauc).
%   D:       asperity size. Function of nominal size per asperity (D0),
%            normal stress (sn), yield stress (sy)
%   tw_ad:   critical weakening time in adiab. conditions. Function of
%            thermal diffusivity (ath), asperity size (D), heat capacity
%            (rhoc), weakening temperature (Tw), initial temperature (T0),
%            contact shear strength (tauc), initial shear stress (tau0),
%            shear zone width (W), number of contacts (Nc)
%   tw_pl:   critical weakening time in s.o.a.p. conditions. Function of
%            thermal diffusivity (ath), asperity size (D), heat capacity
%            (rhoc), weakening temperature (Tw), initial temperature (T0),
%            contact shear strength (tauc), initial shear stress (tau0),
%            number of contacts (Nc)
%   tau_pulse:  shear stress for pulse/crack transition. Function of
%               initial shear strength (tau0), weakening velocity (Vw),
%               shear modulus (mu), shear wave speed (cs)
%   Gc_ad:   max. fracture energy in adiab. conditions. Function of shear
%            zone width (W), heat capacity (rhoc), weakening temperature
%            (Tw), initial temperature (T0)
%   Fp:      function F'(x) appearing in the computation of fracture energy
%            under slip on a plane conditions. Computed numerically and
%            them interpolated.
%   Gc_pl:   fracture energy in s.o.a.p. conditions. Function of asperity
%            size (D), contact shear strength (tauc), initial shear
%            strength (tau0), number of contacts (Nc), normalised time
%            (tnorm).

% parameters and functions for flash heating computations

%% (1) thermal properties of individual minerals

data = dlmread('Thermal_properties_Qz_Olv.txt','\t',1,0);

%quartz, from Clauser and Huenges, 1995
T = data(1:9,1)-273.15;
kT = .5*(data(1:9,2) + data(10:18,2));
rhoc = .5*(data(1:9,4) + data(10:18,4));
ath = .5*(data(1:9,3) + data(10:18,3))/1e6;
%kT = data(10:18,2);
%rhoc = data(10:18,4);
%ath = data(10:18,3)/1e6;

%interpolate data
kT_qz = @(x) interp1(T,kT,x,'spline');
rhoc_qz = @(x) interp1(T,rhoc,x,'spline');
ath_qz = @(x) interp1(T,ath,x,'spline');

%olivine, from Clauser and Huenges, 1995
kT = data(19:27,2);
rhoc = data(19:27,4);
ath = data(19:27,3)/1e6;

kT_olv = @(x) interp1(T,kT,x,'spline');
rhoc_olv = @(x) interp1(T,rhoc,x,'spline');
ath_olv = @(x) interp1(T,ath,x,'spline');

%% (2) thermal properties of rocks, from Vosteen and SchellSchmidt, 2003 

data = dlmread('Thermal_properties_rocks.txt','\t',1,0);
T1 = data(1:10,1);
rhoc1 = data(1:10,2)*1e6;
T2 = data(11:20,1);
rhoc2 = data(11:20,2)*1e6;

%fit rhoc data with a polynomial
%poly_rhoc_rock = polyfit(T2,rhoc2*1e6,1);
%rhoc_rock = @(x) polyval(poly_rhoc_rock,x);

%interpolate data for rhoc
rhoc_rock = @(x) interp1(T2,rhoc2,x,'linear','extrap');
%use their formula for conductivity
kT_rock = @(x) 3.2./(0.99 + x.*(0.003 - 0.0042./3.2));
ath_rock = @(x) kT_rock(x)./rhoc_rock(x);

%% (3) Yield strength data from Hardness

%Quartz data from Evans 1984
data = dlmread('Hardness_Qz_T.txt','\t',3,0);
data([21:26  116:119],:) = [];
T = data(:,1);
sy = data(:,3);

%polynomial fit
poly_sy_qz = polyfit(T,sy*1e9,2);
sy_qz = @(x) polyval(poly_sy_qz,x);

%Olivine data from Evans and Goetze, 1979
data = dlmread('Hardness_Olivine_T.txt','\t',3,0);
T = data(:,1);
sy = data(:,2);

%polynomial fit
poly_sy_olv = polyfit(T,sy*1e9,2);
sy_olv = @(x) polyval(poly_sy_olv,x);

%% flash heating functions

%weakening velocity
Vw = @(ath,D,rhoc,Tw,T0,tauc) pi.*ath./D.*(rhoc.*(Tw-T0)./tauc).^2;
%asperity size
D = @(D0,sn,sy) D0.*sqrt(sn./sy);
%critical time: adiabatic
tw_ad = @(ath,D,rhoc,Tw,T0,tauc,tau0,W,Nc) D.*(W./Nc).*tauc.^2./(tau0.*rhoc.*(Tw-T0).*pi.*ath);
%critical time: plane
tw_pl = @(ath,D,rhoc,Tw,T0,tauc,tau0,Nc) (1./ath).*(D.*tauc.^2./(pi*Nc.*tau0.*rhoc.*(Tw-T0))).^2;
%stress for pulse/crack transition (Zheng and Rice 1998)
tau_pulse = @(tau0,Vw,mu,cs) sqrt(2.*tau0.*Vw.*mu./cs);
%Gc: adiabatic
Gc_ad = @(W,rhoc,Tw,T0) W.*rhoc.*(Tw-T0);
%Gc: slip on a plane
[x,~,Y] = flashheatingsoap_G(15000,2000,300,7);
Fp = @(tn) interp1(x,Y,tn); 
Gc_pl = @(D,tauc,tau0,Nc,tnorm) tauc.^2.*D./(Nc.*pi*tau0) .*Fp(tnorm);