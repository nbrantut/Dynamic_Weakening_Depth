%WATERPROPERTIES
% This routine defines functions to compute water properties. 
% these are computed as finite differences from the
% computed density from the IAPWS standard. Requires the 
% <a href="matlab:web('http://www.peter-junglas.de/fh/water95/index.html')">water95</a> package
% by Peter Junglas, and the viscosity function visc_w(T, rho) also computed
% from IAPWS.
%
% Output (directly in the workspace):
%   lambda_f:  thermal expansivity
%   beta_f:    compressibility
%   eta_f:     viscosity

%thermal expansivity
lambda_f = @(p,T) (-1/density(p,T))*(density(p,T+1)-density(p,T-1))/2;

%compressibility
beta_f = @(p,T) (1/density(p,T))*(density(p+1,T)-density(p-1,T))/2;

%viscosity
eta_f = @(p,T) visc_w(T,density(p,T));