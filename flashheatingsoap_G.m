function [xtot,ytot,Ytot] = flashheatingsoap_G(xmax, Npts,Nskip,imax)
%[xtot,ytot] = flashheatingsoap_G(xmax, Npts,Nskip,imax)
%
%This function computes the fracture energy for flash heating under the
%slip obn a plane approximation. This is not straightforward sinece the
%shear stress history is only known numerically (function
%flashheatingsoap), and a straightforward integral is not accurate at small
%slips. So here the computation is first done on a coarse grid up to a max
%slip given by xmax, for a number of points equal to Npts. Then the first
%Nskip points are skipped, and the G is computed again up to a fraction of
%xmax (given by x(Nskip)). This operation is repeated imax times, with the
%results for the first Nskip points discarded every time and replaced by
%those computed on the finer and finer grid. A good set of values for the
%parameters are: xmamx=15e3; Npts=2000;Nskip=300; imax=7;
%
%input:
%    xmax:    max. nondimensional slip (at V=cte=1m/s) up to which G is
%             computed.
%    Npts:    number of points to make the regular grid in x
%    Nskip:   number of points to discard at small slip (where solution is
%             not accurate).
%    imax:    number of iterations to perform.
%
%ouput:
%    xtot:   nondimensional slip
%    Ytot:   nondimensional fracture energy

[x,y] = flashheatingsoap(xmax,Npts);
Y = cumtrapz(x,(1-y).^2) - x.*(1-y).^2;

xtot = x(Nskip:end);
ytot = y(Nskip:end);
Ytot = Y(Nskip:end);

for i=1:imax
    [x1,y1]=flashheatingsoap(x(Nskip),Npts);
    Y1 = cumtrapz(x1,(1-y1).^2) - x1.*(1-y1).^2;

    
    xtot = [x1(Nskip:end-1); xtot];
    ytot = [y1(Nskip:end-1); ytot];
    Ytot = [Y1(Nskip:end-1); Ytot];
    
    x = x1;
end

%plot(xtot,Ytot)