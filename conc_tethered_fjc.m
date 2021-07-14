function [c,e_cor] = conc_tethered_fjc(L,P,h,s)
%This fucntion is used to calculate the average effective concentration of
% one end of a chain on the membrane. The other end of the chain is anchored at 
% a distance h from the membrane.
 % Inputs:
 %   L: Contour length (nm)
 %   P: Persistence length (nm)
 %   h: Distance of the anchor from the membrane (nm)
 %   s: The area per lipid (typically 0.7 nm*nm).
 % Output:
 %   c: effective concentration on the membrane (M)
 %   e_cor: Energy correction of the binding energy due to the tether (kT)
 
 b=4.*P.*L;  % Calculate <r^2>, since Kuhn length a=2P
 b=3./b;  % b=3/(4PL)
 c=1./s.*sqrt(b./pi).*exp(-h.*h.*b);
 c=c./0.602;  % Convert to M. 0.602=NA*1e-24
 e_cor=-log(c);
end

