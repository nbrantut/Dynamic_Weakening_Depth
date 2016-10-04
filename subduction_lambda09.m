%SUBDUCTION_LAMBDA09
% This routine defines the set of parameters relevant for a subduction zone
% with near lithostatic pore pressure gradient, with fault zone properties
% corresponding to a "granite" fault
% gouge (used as a model for igneous rocks).
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

%slip rate
V=1;

%rate strengthening parameter (a-b)
%from Rice et al 2014, Platt et al 2014
ab = 0.025;

%frictional coefficient
%after flash heating
f = 0.6;

%thermal diffusivity, assumed constant, in m^2/s
ath = 1e-6;

%heat capacity
rhoc = 2.7e6;

%nominal porosity
%n0 = 10e-2;
n0 = @(Peff) 0.05*(1 + exp(-0.022e-6*Peff));

%permeability
%granite gouge after 150 mm slip, Zhang et al JSG 1999
kF = @(Peff) 1e-19*exp(-0.0028*Peff/1e6);

%compressibility of solid grains
%for granite, as quoted in Rice 2006
beta_s = 2e-11;

%drained compressibility
%granite gouge, Zhang et al JSG 1999
%obtained from the dilatation vs Peff fit as
% beta_d = - d (dil./100)/dPeff
%where dil. is the dilatation in % (as plotted)
beta_d = @(Peff) beta_s + 0.022e-6*5.0e-2*exp(-0.022e-6*Peff);

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
zmax = 55;

%pore pressure factor
lambda = 0.9;%32;

%rock density (to compute sigma_n)
rho_rock = 3200;
%gravity
g = 9.8;

%geotherm in C/m
geo = 6e-3;

%% compute arrays from parameter values

%depth in m
depth = linspace(zmin,zmax).*1e3;

%normal stress in Pa
sn = rho_rock*g*depth + 100e6;

%pore pressure in Pa
p0 = 2.8e3*g*depth + 100e6;%1e3*g*depth + 100e6; %2.8e3*g*depth + 100e6;

%temperature in C
T0 = 0 + geo *depth;