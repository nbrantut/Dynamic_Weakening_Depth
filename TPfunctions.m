%TPFUNCTIONS
% This routine defines a set of functions to compute thermal pressurisation
% parameters.
% 
% Output (generated in the workspace):
%   Dc:   slip weakening distance in a.u. conditions. Function of heat 
%         capacity (rhoc), friction coef. (f), thermal 
%         pressurisation parameter (A) and shear zone width (h).
%   gw:   strain weakening distance in a.u. conditions. Function of heat
%         capacity (rhoc), friction coef. (f), thermal pressurisation
%         factor (A).
%   p_au: pore pressure evolution in a.u. conditions. Function of
%         normalised slip (d), normal stress (sn) and initial pore pressure
%         (p0).
%   T_au: temperature evolution in a.u. conditions. Function of normalised
%         slip (d), normal stress (sn), initial pore pressure (p0), initial
%         temperature (T0) and thermal pressurisation factor (A).
%   G_au: fracture energy in a.u. conditons. Function of normalised slip 
%         (d), normal stress (sn), initial pore pressure (p0), friction
%         coef. (f), critical slip weakening distance (dc).
%   Lstar:  critical weakening distance in s.o.a.p. conditions. Function of
%           heat capacity (rhoc), slip rate (V), friction coef. (f), thermal
%           pressurisation factor (A), heat diffusivity (ath), hydraulic
%           diffusivity (ahy).
%   p_soap: pore pressure evolution in s.o.a.p. conditions. Function of
%           normalised slip (d), normal stress (sn) and initial pore
%           pressure (p0).
%   T_soap: temperature evolution in s.o.a.p. conditions. Function of
%           normalised slip (d), normal stress (sn), initial pore pressure
%           (p0), initial temperature (T0), thermal pressurisation factor 
%           (A), thermal diffusivity (ath) and hydraulic diffusivity (ahy).
%   G_soap: fracture energy in s.o.a.p. conditions. Function of normalised
%           slip (d), normal stress (sn), initial pore pressure (p0),
%           friction coefficient (f) and critical distance (lstar).
%   DTmax_au:   max. T rise in a.u. conditions. Function of thermal
%               pressurisation factor (A), normal stress (sn) and initial
%               pore pressure (p0).
%   DTmax_soap: max. T rise in s.o.a.p. conditions. Function of thermal
%               diffusivity (ath), hydraulic diffusivity (ahy), thermal
%               pressurisation factor (A), normal stress (sn), initial pore
%               pressure (p0).
%   Wrsf:   self localising shear zone width (formula from Rice et al.,
%           JGR, 2014; Platt et al., JGR, 2014). Function of heat capacity
%           (rhoc), slip rate (V), thermal pressurisation factor (A),
%           friction coef. (f), thermal diffusivity (ath), hydraulic
%           diffusivity (ahy), rate hardening parameter (ab).


%adiabatic, undrained slip weakening distance
Dc = @(rhoc,f,A,h) rhoc.*h./(A.*f);

%adiabatic undrained weakening strain
gw = @(rhoc,f,A) rhoc./(A.*f);

%pore pressure and temperature solution with normalised slip d = D/Dc
p_au = @(d,sn,p0) p0 + (sn-p0).*(1 - exp(-d));
T_au = @(d,sn,p0,T0,A) T0 + ((sn-p0)/A).*(1 - exp(-d));

%fracture energy with normalised slip d = D/Dc
G_au = @(d,sn,p0,f,dc) f.*(sn-p0).*dc.*(1-(1+d).*exp(-d));

%slip on a plane characteristic weakening distance
Lstar = @(rhoc,V,f,A,ath,ahy) (1./V).*(2.*rhoc.*(sqrt(ath)+sqrt(ahy))./(f.*A)).^2;

%pore pressure and temperature solution with normalised slip d = D/Lstar
p_soap = @(d,sn,p0) sn - (sn-p0).*erfcx(sqrt(d));
T_soap = @(d,sn,p0,T0,A,ath,ahy) T0 + (1+sqrt(ahy/ath)).*((sn-p0)/A).*erfcx(sqrt(d));

%fracture energy with normalised slip d = D/Lstar
G_soap = @(d,sn,p0,f,lstar) f.*(sn-p0).*lstar.*(erfcx(sqrt(d)).*(1-d) -1 + 2*sqrt(d/pi));

%analytical solutions for max temperature
DTmax_au = @(A,sn,p0) (sn-p0)./A;
DTmax_sp = @(ath,ahy,A,sn,p0) (1+sqrt(ahy./ath)) .* (sn-p0)./A;

%strain localisation (l.s.a)
Wrsf = @(rhoc,V,f,A,ath,ahy,ab) pi^2*ab./f.^2.*rhoc./A.*(ahy+ath)./V;