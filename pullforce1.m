function f = pullforce1(D,Lm,Ls,p_m,Emax,dc,d1)
% Calculate the net pulling force between two opposing membranes (the pulling
% force provided by E-Syt transmembrane binding minus the membrane repulsive force)
% and the corresponding total energy

% Inputs:
% D: Membrane separation
% Lm: Contour length corresponding to a C2 binding energy state
% Ls: Hard core size corresponding to Lm. The extension changes required to
% pass the energy barriers are set by this array.
% p_m: Persistence length of the polypeptide
% Emax, dc, d1: parameters for membrane repulsive energy 
% V(d)=Emax*exp[-(d-dc)/d1] with d the membrane distance

% Outputs:
% f: Net membrane pulling force (minus repulsive force) in pN

% All distance parameters have unit of nm, and energy, kT.
     
        % Membrane repulsive energy (kT) and force (pN)
        [~,fm]=membrane_potential_f2(D,Emax,dc,d1);

        r=(D-Ls)./Lm;  % Extension to contour length ratio
        if(r<0)
            r=0;
        end
        [force, ~]=MS(r,p_m);   % Stretching force and energy per contour length
        f=force-fm;  % State pulling force
        
end

