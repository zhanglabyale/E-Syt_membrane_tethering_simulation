function [force,energy] = MS(ratio,P)
% This function calculates the force and entropic energy using the
% Marko-Siggia formula
% ratio:  extension-to-contour length ratio (scalar or vector)
% P:      Persistence length (scalar value), nm
% force:  Stretching force (pN)
% energy: Entropic energy per contour length, nm x pN

% Check if ratio is close to or greater than 1. If it is, set to a pre-set
% maximum value. This value may need change.
rmax=0.95;  
r=ratio;
index=ratio>rmax;
r(index)=rmax;
force=4.1./P.*(0.25./(1-r).^2+r-0.25);
energy=4.1./P./4./(1-r).*r.*r.*(3-2.*r);
force(index)=max(force);
energy(index)=inf;
end

