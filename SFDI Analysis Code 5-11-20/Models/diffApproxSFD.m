%% Standard Diffusion Approximation in the Spatial Frequency Domain
% Calculate normalized diffuse SFD reflectance with the standard
% diffusion approximation.
% adapted from Cuccia et al. 2009
%
function Rd = diffApproxSFD(mua,musp,n,fx)
% diffApproxSFD(mua,musp,n,fx) returns the diffuse reflectance from a 
% homogeneous medium with bulk absorption coefficient (mua), reduced  
% scattering coefficient, refractive index n, and spatial frequency fx.
%
% mua and musp can be vectors or N-D arrays of the same size,
% however, n and fx must be scalar.
%
% mua, musp, and fx have units of [1/mm]
% n is [unitless]
%
% EXAMPLE:
%   Rd = diffApproxSFD([0.001,0.01,0.1],[0.2,0.4,4],1.4,0.1)
%   returns
%   Rd = [0.0310    0.0873    0.4446]

% wavenumber in the x-direction
kx = 2*pi*fx;

% effective reflection coefficient
Reff = 0.0636*n+0.668+(0.71/n)-(1.44/(n^2));

% proportionality constant
A = (1-Reff)/(2*(1+Reff));

% reduced albedo
a = musp./(mua+musp);

% effective interaction coefficient for an SFD source
mueff = ((3*mua.*(mua+musp)) + kx^2).^(1/2);

% numerator
num = 3*A*a;

% denomenator factor 1
den1 = (mueff./(mua+musp)) + 1;

% denomenator factor 2
den2 = (mueff./(mua+musp)) + 3*A;

% final Rd calculation for function return
Rd = num./(den1.*den2);

end