%PROFILE_TP
% This routine computes depth profiles for TP parameters. It does not need
% any input, but the type of profile (continental, oceanic, etc) must be
% entered directly in the file.
%
% This routine makes use of other functions 'waterproperties' and
% 'TPfunctions'.
%
% Output (generated in the workspace):
%   depth:      array of depths (in m)
%   sn:         array of normal stresses (in Pa)
%   p0:         array of pore pressure (in Pa)
%   T0:         array of temperatures (in C)
%   nominal:    a structure containing parameter values computed at 
%               nominal p0, T0 sn conditions, with fields:
%       .kF:          permeability
%       .beta_n_v:    pore space compressibility assuming relaxed off fault
%                     stresses
%       .beta_n_el:   pore space compressibilty assuming elastic off fault
%                     stresses
%       .lambda_n_v:  pore space thermal expansivity assuming relaxed off
%                     fault stresses
%       .lambda_n_el: pore space thermal expansivity assuming elastic off
%                     fault stresses
%       .n0:          porosity
%       .beta_f:      water compressibility
%       .lambda_f:    water thermal expansivity
%       .eta_f:       water viscosity
%       .A:           thermal pressurisation factor (\Lambda)
%       .ahy:         hydraulic diffusivity
%       .Wrsf:        self localised shear zone width
%       .Dc:          adiabatic undrained slip weakening distance (using
%                     W=Wrsf)
%       .Lstar:       slip on a plane slip weakening distance
%       .G_au:        fracture energy (adiabatic undrained, using W=Wrsf)
%       .G_soap:      fracture energy (slip on a plane) at 1 m slip
%       .DTmax_au:    max. T rise (adiabatic undrained)
%       .DTmax_soap:  max. T rise (slip on a plane)
%   pathavgd:   structure containing parameter values computed from
%               path-averaged p, T value, with two fields:
%       .au:          structure obtained from path-averaged values along
%                     the adiabatic, undrained p T path.
%       .soap:        structure obtained from path-averaged values along
%                     the slip on a plane p T path.
%       Both substructures .au and .soap contain the same fields as
%       'nominal', with the exception of the fracture energy. For a.u.
%       conditions, three values of G are computed, using W=10um (G_10um),
%       100um (G_100um) and 1mm (G_1mm). For s.o.a.p conditions, G is
%       computed at 1 m slip only. The max. T rise (au.DTmax and
%       soap.DTmax) is computed using the relevant formulae (a.u. or
%       s.o.a.p.) in each case.
%
%% load parameters:
%this can be either (uncomment the relevant line):
%  continental_granite;
%  continental_claygoue;
%  oceanic;
%  subduction;
%  subduction_lambda09;


%water properties
waterproperties;

%functions to compute thermal pressurisation parameters
TPfunctions;


%% nominal parameters vs depth

%permeability
nominal.kF = kF(sn-p0);

%pore space compressibility
%   assuming relaxed stresses
nominal.beta_n_v  = beta_n_v(sn-p0);
%   assuming elastic stresses
nominal.beta_n_el = beta_n_el(sn-p0);

%pore space thermal expansivity
%   assuming relaxed stresses
nominal.lambda_n_v  = 0*depth + lambda_n_v;
%   assuming elastic stresses
nominal.lambda_n_el = lambda_n_el(sn-p0);

%porosity
nominal.n0 = n0(sn-p0);

%water properties
for k=1:length(depth)
    %compressibility
    nominal.beta_f(k) = beta_f(p0(k),T0(k)+273.15);
    %thermal expansion
    nominal.lambda_f(k) = lambda_f(p0(k),T0(k)+273.15);
    %viscosity
    nominal.eta_f(k) = eta_f(p0(k),T0(k)+273.15);
end

%derived parameters
%TP factor, using elastic approx.
nominal.A = (nominal.lambda_f - nominal.lambda_n_el)./(nominal.beta_f + nominal.beta_n_el);

%hydraulic diffusivity
nominal.ahy = nominal.kF./(nominal.eta_f.*(nominal.beta_f + nominal.beta_n_el).*nominal.n0);

%shear zone width
nominal.Wrsf = Wrsf(rhoc,V,f,nominal.A,ath,nominal.ahy,ab);

%Dc (adiabatic undrained)
nominal.Dc = Dc(rhoc,f,nominal.A,nominal.Wrsf);

%Lstar (slip on a plane)
nominal.Lstar = Lstar(rhoc,V,f,nominal.A,ath,nominal.ahy);

%Fracture energy for 1m slip
%   under adiabatic, undrained conditions
nominal.G_au = G_au(1./nominal.Dc,sn,p0,f,nominal.Dc);
%   under slip on plance conditions
nominal.G_soap = G_soap(1./nominal.Lstar,sn,p0,f,nominal.Lstar);

%Max temperature rise
%   under adiabatic undrained conditions
nominal.DTmax_au = DTmax_au(nominal.A,sn,p0);
%   under slip on a plane conditions
nominal.DTmax_soap = DTmax_sp(ath,nominal.ahy,nominal.A,sn,p0);

%% path averaged parameters with depth

%the approach to get path averaged parameters is the following:
%(1) compute p and T evolution during slip using nominal values of
%parameters at each depth, and for each scenario (a.u. or s.o.a.p.). The
%slip distance over which we compute is from 0 to either Dc (a.u.) or Lstar
%(s.o.a.p.), that is, a normalised slip from 0 to 1. We choose this slip
%interval because this is the relevant one at the onset of slip (the one
%that is used to define a slip weakening distance).
%(2) then compute averaged parameter values along that path, according to
%
%   X_m = \int_{0}^{1} X[p(s),T(s)] ds
%
%where X_m is the averaged value of parameter X, and s is slip.

%normalised slip
d = linspace(0,1);
N = length(d);

for k=1:length(depth)
    
    %compute p, T paths:
    %   adiabatic, undrained
    p_path_au = p_au(d,sn(k),p0(k));
    T_path_au = T_au(d,sn(k),p0(k),T0(k),nominal.A(k));
    %   slip on a plane
    p_path_soap = p_soap(d,sn(k),p0(k));
    T_path_soap = T_soap(d,sn(k),p0(k),T0(k),nominal.A(k),ath,nominal.ahy(k));

    %permeability
    pathavgd.au.kF(k) = (1/N)*trapz(kF(sn(k) - p_path_au));
    pathavgd.soap.kF(k) = (1/N)*trapz(kF(sn(k) - p_path_soap));
    
    %porosity
    pathavgd.au.n0(k) = (1/N)*trapz(n0(sn(k) - p_path_au));
    pathavgd.soap.n0(k) = (1/N)*trapz(n0(sn(k) - p_path_soap));
    
    %pore space compressibility
    pathavgd.au.beta_n_v(k) = (1/N)*trapz(beta_n_v(sn(k) - p_path_au));
    pathavgd.soap.beta_n_v(k) = (1/N)*trapz(beta_n_v(sn(k) - p_path_soap));
    
    pathavgd.au.beta_n_el(k) = (1/N)*trapz(beta_n_el(sn(k) - p_path_au));
    pathavgd.soap.beta_n_el(k) = (1/N)*trapz(beta_n_el(sn(k) - p_path_soap));
    
    %pore space thermal expansivity
    pathavgd.au.lambda_n_v(k) = lambda_n_v;
    pathavgd.soap.lambda_n_v(k) = lambda_n_v;
    
    pathavgd.au.lambda_n_el(k) = (1/N)*trapz(lambda_n_el(sn(k)-p_path_au));
    pathavgd.soap.lambda_n_el(k) = (1/N)*trapz(lambda_n_el(sn(k)-p_path_soap));
    
    %water properties:
    %initialise values at each slip for both scenarios
    au_tbeta_f = 0*d;
    au_tlambda_f = 0*d;
    au_teta_f = 0*d;
    
    soap_tbeta_f = 0*d;
    soap_tlambda_f = 0*d;
    soap_teta_f = 0*d;
    
    %compute water properties at each displacement
    for l=1:N
        %adiabatic, undrained scenario
        au_tbeta_f(l) = beta_f(p_path_au(l),T_path_au(l)+273.15);
        au_tlambda_f(l) = lambda_f(p_path_au(l),T_path_au(l)+273.15);
        au_teta_f(l) = eta_f(p_path_au(l),T_path_au(l)+273.15);
        %slip on a plane scenario
        soap_tbeta_f(l) = beta_f(p_path_soap(l),T_path_soap(l)+273.15);
        soap_tlambda_f(l) = lambda_f(p_path_soap(l),T_path_soap(l)+273.15);
        soap_teta_f(l) = eta_f(p_path_soap(l),T_path_soap(l)+273.15);
    end
    
    %now compute averages
    pathavgd.au.beta_f(k) = (1/N)*trapz(au_tbeta_f);
    pathavgd.au.lambda_f(k) = (1/N)*trapz(au_tlambda_f);
    pathavgd.au.eta_f(k) = (1/N)*trapz(au_teta_f);
    
    pathavgd.soap.beta_f(k) = (1/N)*trapz(soap_tbeta_f);
    pathavgd.soap.lambda_f(k) = (1/N)*trapz(soap_tlambda_f);
    pathavgd.soap.eta_f(k) = (1/N)*trapz(soap_teta_f);
end
    
%% derived parameters
%TP factor, using elastic approx.
pathavgd.au.A = (pathavgd.au.lambda_f - pathavgd.au.lambda_n_el)./(pathavgd.au.beta_f + pathavgd.au.beta_n_el);
pathavgd.soap.A = (pathavgd.soap.lambda_f - pathavgd.soap.lambda_n_el)./(pathavgd.soap.beta_f + pathavgd.soap.beta_n_el);

%hydraulic diffusivity
pathavgd.au.ahy = pathavgd.au.kF./(pathavgd.au.eta_f.*(pathavgd.au.beta_f + pathavgd.au.beta_n_el).*pathavgd.au.n0);
pathavgd.soap.ahy = pathavgd.soap.kF./(pathavgd.soap.eta_f.*(pathavgd.soap.beta_f + pathavgd.soap.beta_n_el).*pathavgd.soap.n0);

%shear zone width
pathavgd.au.Wrsf = Wrsf(rhoc,V,f,pathavgd.au.A,ath,pathavgd.au.ahy,ab);
pathavgd.soap.Wrsf = Wrsf(rhoc,V,f,pathavgd.soap.A,ath,pathavgd.soap.ahy,ab);

%undrained adiabatic
pathavgd.au.gw = gw(rhoc,f,pathavgd.au.A);

%Lstar
pathavgd.soap.Lstar = Lstar(rhoc,V,f,pathavgd.soap.A,ath,pathavgd.soap.ahy);

%Gc for 1m slip
pathavgd.au.G_10um = G_au(1./(1e-5.*pathavgd.au.gw),sn,p0,f,1e-5.*pathavgd.au.gw);
pathavgd.au.G_100um = G_au(1./(1e-4.*pathavgd.au.gw),sn,p0,f,1e-4.*pathavgd.au.gw);
pathavgd.au.G_1mm = G_au(1./(1e-3.*pathavgd.au.gw),sn,p0,f,1e-3.*pathavgd.au.gw);
pathavgd.soap.G = G_soap(1./pathavgd.soap.Lstar,sn,p0,f,pathavgd.soap.Lstar);

%DTmax
pathavgd.au.DTmax = DTmax_au(pathavgd.au.A,sn,p0);
pathavgd.soap.DTmax = DTmax_sp(ath,pathavgd.soap.ahy,pathavgd.soap.A,sn,p0);
