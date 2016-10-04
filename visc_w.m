function mu = visc_w(T,rho)
%mu = VISC_W(T,rho);
%
%calculates the viscosity of water (mu) as a function of temperature and
%density (SI units).
%
%computed from <a href="matlab:web('http://iapws.org/relguide/visc.pdf')">IAPWS 2008 release</a>.
%
%CRITICAL ENHANCEMENT NOT INCLUDED (result not accurate near critical
%point).
%
%input:
%   T:      temperature in K
%   rho:    density is kg/m3
%output:
%   mu:     viscosity in Pa s
%
%code written by Nicolas Brantut, University College London, Seismolab/RIPL
%last modified: 15 Apr. 2014

%% reference values
T0 = 647.096;
rho0 = 322.0;
p0 = 22.064;
mu0 = 1.00e-6;

%% constants
H0 = 1.67752;
H1 = 2.20462;
H2 = 0.6366564;
H3 = -0.241605;

H00 = 5.20094e-1;
H10 = 8.50895e-2;
H20 = -1.08374;
H30 = -2.89555e-1;
H01 = 2.22531e-1;
H11 = 9.99115e-1;
H21 = 1.88797;
H31 = 1.26613;
H51 = 1.20573e-1;
H02 = -2.81378e-1;
H12 = -9.06851e-1;
H22 = -7.72479e-1;
H32 = -4.89837e-1;
H42 = -2.57040e-1;
H03 = 1.61913e-1;
H13 = 2.57399e-1;
H04 = -3.25372e-2;
H34 = 6.98452e-2;
H45 = 8.72102e-3;
H36 = -4.35673e-3;
H56 = -5.93264e-4;

%% normalise inputs
t = T/T0;
r = rho/rho0;

%% calculate viscosity factors
m0 = 100*sqrt(t)./(H0 + H1./t + H2./t.^2 + H3./t.^3);

m1 = exp(r.*(1*(H00 + H01.*(r-1) + H02.*(r-1).^2 + H03.*(r-1).^3 + H04.*(r-1).^4) + ...
    (1./t-1).^1.*(H10 + H11.*(r-1) + H12.*(r-1).^2 + H13.*(r-1).^3) + ...
    (1./t-1).^2.*(H20 + H21.*(r-1) + H22.*(r-1).^2) + ...
    (1./t-1).^3.*(H30 + H31.*(r-1) + H32.*(r-1).^2 + H34.*(r-1).^4 + H36.*(r-1).^6) + ...
    (1./t-1).^4.*(H42.*(r-1).^2 + H45.*(r-1).^5) + ...
    (1./t-1).^5.*(H51.*(r-1) + H56.*(r-1).^6) ));

m2 = 1; %no critical point correction

%% calculate viscosity

mu = mu0.*m0.*m1.*m2;