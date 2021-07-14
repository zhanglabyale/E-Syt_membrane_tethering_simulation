
function [V,f]=membrane_potential_f2(d,Emax,dc,d1)
% The opposing energy of bringing two membranes at distance d
% d: distance, vector or scalar; in unit of nm
% Emax: Maximum energy potential at d=dc; unit is kT
% dc: starting point for expotential decay above this point
% d1: the decaying length of the membrane repulsive potential. Unit nm
% V: membrane potential in unit kT
% potentia V(d)=Emax*exp(-(d-dc)/d1)
% force f=Emax*exp(-(d-dc)/d1)/d1 in the unit of pN

V=zeros(size(d));
f=zeros(size(d));
kT=4.1;
ds=d-dc;
V=Emax.*exp(-ds./d1);
f=kT.*V./d1;

end


