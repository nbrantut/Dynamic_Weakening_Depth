function [x,y] = flashheatingsoap(x_max,N)
%[x,y] = flashheatingsoap(x_max,N)
%
%this function computes an numerical solution to the integral equation for
%temperature evolution during weakening by flash heating in the slip on a
%plane case. It is all non dimensionalised so the equation that is solved
%is the following:
%
%   y(x) = 1/2\sqrt(pi) \int_0^x (1- y(t))^2/sqrt(x-t) dt
%
%which is equivalent to (after integration by parts to regularise the numerical problem):
%
%   0 = -y(x) + sqrt(x/pi) [1 - y(x)]^2 + 1/2\sqrt(pi) \int_0^x ...
%   ... (2y(x) - y(x)^2 - 2y(t) + y(t)^2)/sqrt(x-t) dt
%
%For small x, the solution for y is asymptotically approximated by
%
% y(x) = (1/2)(1 - exp(x)erfc(\sqrt(x)))
%
%Input:
%   x_max:  the maximum x up to which the function is computed
%   N:      the number of points between 0 and x_max, regularly spaced
%
%Output:
%   x:      the array of x's: x=linspace(0,x_max,N+1);
%   y:      the value of y(x).


%initial guess for y
y0 = zeros(N+1,1);
%step and x array
h = x_max/N;
x = (0:h:x_max)';

%grids to compute matrices efficiently
xxj = meshgrid(x);
xxi = xxj';
sqxx = tril(1./sqrt(xxi - xxj),-1);

%arrays with indices to compute matrices efficiently
ind = 1:N+1;
indj = meshgrid(ind);
indi = indj';

%weight function for trapezoidal rule to approximate the integral
w = ones(N+1);
w(:,1) = 1/2;
w = tril(w,-1);

%thats i a useful constant
const = h/(2*sqrt(pi)).*sum(w.*sqxx,2);

%initialise
y=y0;

%set arbitrary large dy to start the loop
dy=1;

%that's the tolerance
tol = 1e-6;

while max(abs(dy))>tol
    
    %residual
    f = -y + sqrt(x/pi).*(1-y).^2 + h/(2*sqrt(pi)).*sum(w.*sqxx.*(2.*y(indi) - y(indi).^2 - 2.*y(indj) + y(indj).^2),2);

    %jacobian
    J = diag(-1 -2.*sqrt(x/pi).*(1-y) + 2.*(1-y).*const) - h/sqrt(pi).*w.*sqxx.*(1-y(indj));

    %these two f and J are for the linear case, for which the small x
    %asymptotic solution is valid. It was used as a test of the numerical
    %method for benchmarking.
    %f = -y +sqrt(x/pi).*(1-2.*y) + h/sqrt(pi).*sum(w.*sqxx.*(y(indi) - y(indj)),2);
    %J = diag(-1 -2*sqrt(x/pi)+ 2.*const) - h/sqrt(pi).*w.*sqxx;

    %compute Newton raphson step
    dy = -J\f;

    %update y
    y = y + dy;
    
end