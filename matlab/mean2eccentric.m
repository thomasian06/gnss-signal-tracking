function [E] = mean2eccentric(M,ecc) 

%==========================================================================
%==========================================================================
% [E] = mean2eccentric(M,e)
%
% Calculates the eccentric anomaly from the mean anomaly using Newton's 
%  method. The tolerance is set to 1e-12.
%
%
% Author: Ben K. Bradley
% Last Revision Date: 26 October 2010
%
%
%
% INPUT:            Description                                       Units
%
%  M          - mean anomaly (single value or vector)                   rad
%  ecc        - eccentricity orbit (single value or vector)      
%
% OUTPUT:
%
%  E          - eccentric anomaly (same size as input M)                rad
%
%
% Coupling:
%
%  none
%
% References:
%
%  [1] Curtis, H.D. "Orbital Mechanics for Engineering Students".
%          Elsevier Ltd. 2005.
%
%==========================================================================
%==========================================================================

if any(find((ecc > 0.999))) || any(find((ecc < 0)))
    errordlg({'Eccentric anomaly not solved for in mean2eccentric.m','Not an elliptic orbit.'},'Error!');
end

%==========================================================================

% Set tolerance
twopi   = 2*pi;
tol     = 1.0e-12;
maxiter = 20;


% Make sure mean anomaly is positive and within 2pi
M = rem(M,twopi);

M = (M<0)*twopi + M;


%  Make an initial guess for E, Eccentric Anomaly
% if (M < pi)
%     E = M + ecc/2;
% else
%     E = M - ecc/2;
% end
sinM = sin(M);

E = M + (ecc.*sinM) ./ (1 - sin(M + ecc) + sinM);


% Initialize iteration count and error
iter = 1;
err  = 1;


% Iterate to find E, eccentric anomaly
while any(find(abs(err) > tol)) && (iter <= maxiter)
    
    err  = (E - ecc.*sin(E) - M) ./ (1 - ecc.*cos(E));

    E    = E - err;
    
    iter = iter + 1;
    
    if (iter > maxiter)
        warning('Iterations maxed out in mean2eccentric.m'); %#ok<WNTAG>
    end
    
end





