%CONTINENTAL_CLAYGOUGE
% This routine defines the set of parameters relevant for a continental
% crust, with fault zone properties corresponding to a clay-rich fault gouge.
%
% ALL NUMBERS IN SI UNITS
%
% Output (generated in the workspace):
%   V:      slip rate.
%   ab:     rate strengthening parameter
%   f:      friction coefficient
%   ath:    heat diffusivity
%   rhoc:   heat capacity
%   n0:     porosity (function of effective pressure)
%   kF:     permeability (function of effective pressure)
%   beta_s: compressibility of solid
%   beta_d: drained compressibility (function of effective pressure)
%   beta_n_v:    pore space compressibility assuming relaxed off fault
%                stresses (function of effective pressure)
%   beta_n_el:   pore space compressibilty assuming elastic off fault
%                stresses (function of effective pressure)
%   lambda_n_v:  pore space thermal expansivity assuming relaxed off
%                fault stresses (function of effective pressure)
%   lambda_n_el: pore space thermal expansivity assuming elastic off
%                fault stresses (function of effective pressure)
%   zmin, zmax:  depth range
%   lambda:      pore pressure factor
%   rho_rock:    crustal rock density
%   g:           gravity acceleration
%   qT:          surface heat flux
%   kT:          crustal rock conductivity
%   A0:          nominal radioactive heat production
%   hr:          characteristic decay depth for radioactive heat production
%   depth:   array of depths
%   sn:      array of normal stresses
%   p0:      array of initial pore pressures
%   T0:      array of initial temperatures

%% parameters and functions to compute TP, FH

%slip rate (m/s)
V=1;

%rate strengthening parameter (a-b)
%from Rice et al 2014, Platt et al 2014
ab = 0.025;

%frictional coefficient
%after flash heating
f = 0.6;

%thermal diffusivity, assumed constant, in m^2/s
ath = 1e-6;

%heat capacity (MPa/C)
rhoc = 2.7e6;

%nominal porosity (an average representative value from the ones quoted in
%Rice 2006)
n0 = @(Peff) 6.0e-2*exp(-0.0038e-6.*Peff);

%permeability
%MTL gouge, as in Rice 2006
kF = @(Peff) 2.12e-19*exp(-0.0288*Peff/1e6);
%MTL gouge, as in Wibberley and Shimamoto, 2005
%kF = @(Peff) 8.71e-21 * exp(-0.0326*Peff/1e6);

%compressibility of solid grains
%for gouge, as quoted in Rice 2006
beta_s = 1.6e-11;

%drained compressibility
%MTL gouge, as in Rice 2006
beta_d = @(Peff) beta_s + 1.39e-10*exp(-0.00691e-6*Peff);


%thermal expansion of solid grains
%for granite, as quoted in Rice 2006
lambda_s = 2.4e-5;

%compressibility and thermal exp of pore space, assuming relaxed stresses
% eq. A2 in Rice 2006
beta_n_v = @(Peff) (beta_d(Peff) - beta_s)/n0(Peff);
lambda_n_v = lambda_s;

% assuming unrelaxed elastic stresses
% eq. A8 in Rice 2006, using r=1
beta_n_el = @(Peff) (beta_d(Peff) - beta_s).*(beta_d(Peff)+beta_s)./(2*n0(Peff).*beta_d(Peff)) - beta_s;
lambda_n_el = @(Peff) lambda_s.*(1-.5*(beta_d(Peff) - beta_s)./(n0(Peff).*beta_d(Peff)));

% parameters to compute profiles in stress, pore pressure, temperature

%depth interval in km
zmin = 2;
zmax = 25;

%pore pressure factor
lambda = 0.32;

%rock density (to compute sigma_n)
rho_rock = 2800;
%gravity
g = 9.8;

%surface heat flux in W/m^2
%source: Lachenbruch and Sass, 1980, for the SAF
qT = 80e-3; 
%rock average thermal conductivity in W/K/m
%source: Jaupart and Mareschal, 2007; Chapman 1986.
kT = 3;

%% compute arrays from parameter values

%depth in m
depth = linspace(zmin,zmax).*1e3;

%normal stress in Pa
sn = rho_rock*g*depth;

%pore pressure in Pa
p0 = 1e3*g*depth;

%geotherm:
A0=2e-6;
hr=10e3;
T0 = 13 + (A0*hr^2/kT)*(1-exp(-depth/hr)) + (qT-A0*hr)*depth/kT;