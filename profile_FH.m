%PROFILE_TP
% This routine computes depth profiles for FH parameters. It does not need
% any input, but the type of profile (continental, oceanic, etc) must be
% entered directly in the file.
%
% This routine makes use of other function 'FHfunctions'. The computations 
% are stored into an array of structures called "fh", with hopefully
% self-explanatory fields:
%   type:       the geological setting (string)
%   T0:         ambient temperature profile
%   sn:         ambient normal stress profile
%   p0:         ambient pore pressure profile
%   depth:      array with depths
%   D:          asperity diameter
%   Tw:         weakening temperature, assumed =1000C always
%   Vw:         weakening velocity
%   d:          particle (contact) spacing in the gouge
%   tw_ad:      weakening time for adiabatic case
%   tw_pl:      ...................slip on a plane
%   tau_pulse:  min. stress for slip pulse.
%   Gc_ad_10um: Gc for adibatic conditions, W=10um
%   Gc_ad_100um:..............................100um
%   Gc_ad_1mm:  ..............................1mm
%   Gc_pl:      Gc for slip on a plane conditions, 1m slip.

%ALL IN SI UNITS

%% load fuctions and parameters
FHfunctions;

%% continental granite

i=1;
%load geotherms
continental_granite;

fh(i).type = 'continental_granite';
%temperature
fh(i).T0 = T0;
%normal stress (from lithostatic)
fh(i).sn = sn;
%pore pressure (from hydrostatic)
fh(i).p0 = p0;
%depth
fh(i).depth = depth;
%asperity size
fh(i).D = D(200e-6,sn-p0,sy_qz(T0));
%weakening T
fh(i).Tw = 1000;
%gouge width
fh(i).W = 100e-6;
%Number of contacts across the gouge width (needed for tw calculation)
fh(i).Nc = 10;
%weakening velocity
fh(i).Vw = fh(i).Nc.*Vw(ath_qz(T0),fh(i).D,rhoc_qz(T0),fh(i).Tw,T0,0.6*sy_qz(T0));
%for adiabatic model for gouge:
fh(i).tw_ad = tw_ad(ath_qz(T0),fh(i).D,rhoc_qz(T0),fh(i).Tw,T0,0.6*sy_qz(T0),0.6*(sn-p0),fh(i).W,fh(i).Nc);
%for slip on a plane model for bare rocks:
fh(i).tw_pl = tw_pl(ath_qz(T0),fh(i).D,rhoc_qz(T0),fh(i).Tw,T0,0.6*sy_qz(T0),0.6*(sn-p0),fh(i).Nc);
%minimum stress to allow crack-like propagation
fh(i).tau_pulse = tau_pulse(0.6*(sn-p0),fh(i).Vw,30e9,3e3);
%gc adiab
fh(i).Gc_ad_10um = Gc_ad(10e-6,rhoc_qz(T0),fh(i).Tw,fh(i).T0);
fh(i).Gc_ad_100um = Gc_ad(100e-6,rhoc_qz(T0),fh(i).Tw,fh(i).T0);
fh(i).Gc_ad_1mm = Gc_ad(1e-3,rhoc_qz(T0),fh(i).Tw,fh(i).T0);
%gc slip on a plane (slip = 1m)
fh(i).Gc_pl = Gc_pl(fh(i).D,0.6*sy_qz(T0),0.6*(sn-p0),fh(i).Nc,1./fh(i).tw_pl);

%plot
figure;
subplot 221
plot(fh(i).D*1e6, -fh(i).depth/1e3);
xlabel('asperity diameter, {\itD}(\mum)')
ylabel('depth (km)')
title('Continental crust, Quartz')

subplot 222
plot(fh(i).Vw, -fh(i).depth/1e3);
xlabel('weakening velocity, {\itV}_w (m/s)')
ylabel('depth (km)')
title('Continental crust, Quartz')

subplot 223
semilogx(fh(i).tw_ad, -fh(i).depth/1e3,fh(i).tw_pl, -fh(i).depth/1e3);
legend('adiabatic','on plane','location','SouthEast');
xlabel('weakening time, {\itt}_w^{ad} or {\itt}_w^{pl} (s)')
ylabel('depth (km)')
title('Continental crust, Quartz')

subplot 224
plot(fh(i).tau_pulse/1e6, -fh(i).depth/1e3,0.6*(fh(i).sn-fh(i).p0)/1e6,-fh(i).depth/1e3);
legend('\tau_{pulse}','\tau_0','location','NorthEast');
xlabel('stress, \tau_0 and \tau_{pulse} (MPa)')
ylabel('depth (km)')
title('Continental crust, Quartz')

%% oceanic

i=2;
%load geotherms
oceanic;

fh(i).type = 'oceanic';
fh(i).T0 = T0;
fh(i).sn = sn;
fh(i).p0 = p0;
fh(i).depth = depth;
fh(i).D = D(200e-6,sn-p0,sy_olv(T0));
fh(i).Tw = 1000;
%gouge width
fh(i).W = 100e-6;
%Number of contacts across the gouge width (needed for tw calculation)
fh(i).Nc = 10;
%weakening velocity
fh(i).Vw = fh(i).Nc.*Vw(ath_olv(T0),fh(i).D,rhoc_olv(T0),fh(i).Tw,T0,0.6*sy_olv(T0));
%for adiabatic model for gouge:
fh(i).tw_ad = tw_ad(ath_olv(T0),fh(i).D,rhoc_olv(T0),fh(i).Tw,T0,0.6*sy_olv(T0),0.6*(sn-p0),fh(i).W,fh(i).Nc);
%for slip on a plane model for bare rocks:
fh(i).tw_pl = tw_pl(ath_olv(T0),fh(i).D,rhoc_olv(T0),fh(i).Tw,T0,0.6*sy_olv(T0),0.6*(sn-p0),fh(i).Nc);
%minimum stress to allow crack-like propagation
fh(i).tau_pulse = tau_pulse(0.6*(sn-p0),fh(i).Vw,30e9,3e3);
%gc adiab
fh(i).Gc_ad_10um = Gc_ad(10e-6,rhoc_olv(T0),fh(i).Tw,fh(i).T0);
fh(i).Gc_ad_100um = Gc_ad(100e-6,rhoc_olv(T0),fh(i).Tw,fh(i).T0);
fh(i).Gc_ad_1mm = Gc_ad(1e-3,rhoc_olv(T0),fh(i).Tw,fh(i).T0);
%gc slip on a plane (slip = 1m)
fh(i).Gc_pl = Gc_pl(fh(i).D,0.6*sy_olv(T0),0.6*(sn-p0),fh(i).Nc,1./fh(i).tw_pl);

%plot
figure;
subplot 221
plot(fh(i).D*1e6, -fh(i).depth/1e3);
xlabel('asperity diameter, {\itD}(\mum)')
ylabel('depth (km)')
title('Oceanic crust, Olivine')

subplot 222
plot(fh(i).Vw, -fh(i).depth/1e3);
xlabel('weakening velocity, {\itV}_w (m/s)')
ylabel('depth (km)')
title('Oceanic crust, Olivine')

subplot 223
semilogx(fh(i).tw_ad, -fh(i).depth/1e3,fh(i).tw_pl, -fh(i).depth/1e3);
legend('adiabatic','on plane','location','SouthEast');
xlabel('weakening time, {\itt}_w^{ad} or {\itt}_w^{pl} (s)')
ylabel('depth (km)')
title('Oceanic crust, Olivine')

subplot 224
plot(fh(i).tau_pulse/1e6, -fh(i).depth/1e3,0.6*(fh(i).sn-fh(i).p0)/1e6,-fh(i).depth/1e3);
legend('\tau_{pulse}','\tau_0','location','NorthEast');
xlabel('stress, \tau_0 and \tau_{pulse} (MPa)')
ylabel('depth (km)')
title('Oceanic crust, Olivine')

%% subduction

i=3;
%load geotherms
subduction;

fh(i).type = 'subduction';
fh(i).T0 = T0;
fh(i).sn = sn;
fh(i).p0 = p0;
fh(i).depth = depth;
fh(i).D = D(200e-6,sn-p0,sy_olv(T0));
fh(i).Tw = 1000;
%gouge width
fh(i).W = 100e-6;
%Number of contacts across the gouge width (needed for tw calculation)
fh(i).Nc = 10;
%weakening velocity
fh(i).Vw = fh(i).Nc.*Vw(ath_olv(T0),fh(i).D,rhoc_olv(T0),fh(i).Tw,T0,0.6*sy_olv(T0));
%for adiabatic model for gouge:
fh(i).tw_ad = tw_ad(ath_olv(T0),fh(i).D,rhoc_olv(T0),fh(i).Tw,T0,0.6*sy_olv(T0),0.6*(sn-p0),fh(i).W,fh(i).Nc);
%for slip on a plane model for bare rocks:
fh(i).tw_pl = tw_pl(ath_olv(T0),fh(i).D,rhoc_olv(T0),fh(i).Tw,T0,0.6*sy_olv(T0),0.6*(sn-p0),fh(i).Nc);
%minimum stress to allow crack-like propagation
fh(i).tau_pulse = tau_pulse(0.6*(sn-p0),fh(i).Vw,30e9,3e3);
%gc adiab
fh(i).Gc_ad_10um = Gc_ad(10e-6,rhoc_olv(T0),fh(i).Tw,fh(i).T0);
fh(i).Gc_ad_100um = Gc_ad(100e-6,rhoc_olv(T0),fh(i).Tw,fh(i).T0);
fh(i).Gc_ad_1mm = Gc_ad(1e-3,rhoc_olv(T0),fh(i).Tw,fh(i).T0);
%gc slip on a plane (slip = 1m)
fh(i).Gc_pl = Gc_pl(fh(i).D,0.6*sy_olv(T0),0.6*(sn-p0),fh(i).Nc,1./fh(i).tw_pl);

%plot
figure;
subplot 221
plot(fh(i).D*1e6, -fh(i).depth/1e3);
xlabel('asperity diameter, {\itD}(\mum)')
ylabel('depth (km)')
title('Subduction zone, Olivine')

subplot 222
plot(fh(i).Vw, -fh(i).depth/1e3);
xlabel('weakening velocity, {\itV}_w (m/s)')
ylabel('depth (km)')
title('Subduction zone, Olivine')

subplot 223
semilogx(fh(i).tw_ad, -fh(i).depth/1e3,fh(i).tw_pl, -fh(i).depth/1e3);
legend('adiabatic','on plane','location','SouthEast');
xlabel('weakening time, {\itt}_w^{ad} or {\itt}_w^{pl} (s)')
ylabel('depth (km)')
title('Subduction zone, Olivine')

subplot 224
plot(fh(i).tau_pulse/1e6, -fh(i).depth/1e3,0.6*(fh(i).sn-fh(i).p0)/1e6,-fh(i).depth/1e3);
legend('\tau_{pulse}','\tau_0','location','NorthEast');
xlabel('stress, \tau_0 and \tau_{pulse} (MPa)')
ylabel('depth (km)')
title('Subduction zone, Olivine')


%% subduction, near lithostatic pore pressure

i=4;
%load geotherms
subduction_lambda09;

fh(i).type = 'subduction_lambda=0.9';
fh(i).T0 = T0;
fh(i).sn = sn;
fh(i).p0 = p0;
fh(i).depth = depth;
fh(i).D = D(100e-6,sn-p0,sy_olv(T0));
fh(i).Tw = 1000;
%gouge width
fh(i).W = 100e-6;
%Number of contacts across the gouge width (needed for tw calculation)
fh(i).Nc = 10;
%weakening velocity for gouge
fh(i).Vw = fh(i).Nc.*Vw(ath_olv(T0),fh(i).D,rhoc_olv(T0),fh(i).Tw,T0,0.6*sy_olv(T0));
%for adiabatic model for gouge:
fh(i).tw_ad = tw_ad(ath_olv(T0),fh(i).D,rhoc_olv(T0),fh(i).Tw,T0,0.6*sy_olv(T0),0.6*(sn-p0),fh(i).W,fh(i).Nc);
%for slip on a plane model for bare rocks:
fh(i).tw_pl = tw_pl(ath_olv(T0),fh(i).D,rhoc_olv(T0),fh(i).Tw,T0,0.6*sy_olv(T0),0.6*(sn-p0),fh(i).Nc);
%minimum stress to allow crack-like propagation
fh(i).tau_pulse = tau_pulse(0.6*(sn-p0),fh(i).Vw,30e9,3e3);
%gc adiab
fh(i).Gc_ad_10um = Gc_ad(10e-6,rhoc_olv(T0),fh(i).Tw,fh(i).T0);
fh(i).Gc_ad_100um = Gc_ad(100e-6,rhoc_olv(T0),fh(i).Tw,fh(i).T0);
fh(i).Gc_ad_1mm = Gc_ad(1e-3,rhoc_olv(T0),fh(i).Tw,fh(i).T0);
%gc slip on a plane (slip = 1m)
fh(i).Gc_pl = Gc_pl(fh(i).D,0.6*sy_olv(T0),0.6*(sn-p0),fh(i).Nc,1./fh(i).tw_pl);

%plot
figure;
subplot 221
plot(fh(i).D*1e6, -fh(i).depth/1e3);
xlabel('asperity diameter, {\itD}(\mum)')
ylabel('depth (km)')
title('Subduction zone, Olivine')

subplot 222
plot(fh(i).Vw, -fh(i).depth/1e3);
xlabel('weakening velocity, {\itV}_w (m/s)')
ylabel('depth (km)')
title('Subduction zone, Olivine')

subplot 223
semilogx(fh(i).tw_ad, -fh(i).depth/1e3,fh(i).tw_pl, -fh(i).depth/1e3);
legend('adiabatic','on plane','location','SouthEast');
xlabel('weakening time, {\itt}_w^{ad} or {\itt}_w^{pl} (s)')
ylabel('depth (km)')
title('Subduction zone, Olivine')

subplot 224
plot(fh(i).tau_pulse/1e6, -fh(i).depth/1e3,0.6*(fh(i).sn-fh(i).p0)/1e6,-fh(i).depth/1e3);
legend('\tau_{pulse}','\tau_0','location','NorthEast');
xlabel('stress, \tau_0 and \tau_{pulse} (MPa)')
ylabel('depth (km)')
title('Subduction zone, Olivine')