function [fp,fm,e] = force_energy2(D,kmax,kon,koff,Lm,Ls,p_m,s,Emax,dc,d1)
% Calculate the state-dependent E-Syt tethering force and membrane repulsive force, and the
% total energy, given a membrane distance and binding parameters of C2
% domains. Note that here kon, koff, Lm, Ls are one-dimensional array for a single
% E-Syt molecule

% Inputs:
% D: Membrane separation
% kmax: Maximum transition rate at zero energy barrier
% kon, koff: Membrane binding and unbinding rates of two C2 modules (assending order of
% corresponding membrane tether
% Lm: Contour length corresponding to C2 binding energy landscape. The
% unbound state is set to the maximum of the contour lengths for the stable
% bound state
% Ls: Hard core size corresponding to Lm. The extension changes required to
% pass the energy barriers are set by this array.
% p_m: Persistence length of the polypeptide
% s: Area per lipid in the membrane.
% Emax, dc, d1: parameters for membrane repulsive energy. See
% function membrane_potential_f2.m

% Outputs:
% fp: State-dependent protein pulling force in pN (excluding transition
% states)
% fm: D-dependent membrane repulsive force
% e: State-dependent total energy (Memebrane binding, linker stretching, and membrane repulsion) 

% All distance parameters have unit of nm, and energy, kT.
     
        % Membrane repulsive energy (kT) and force (pN)
        [V,fm]=membrane_potential_f2(D,Emax,dc,d1);

        % Calculate the binding energy of the ternminal C2 binding module
        % at zero force or zero membrane separation
        [c,~] = conc_tethered_fjc(Lm(3),p_m,0,s);
        V_one=log(kon(2).*c/koff(2));  % Unbinding energy of the C2E in kT at minimum membrane separation
        
        % Calculate the binding energy of the first C2 binding module at a zero force given
        % that the terminal C2 module is bound.
        L_two=Lm(3)-Lm(1);
        [c,~] = conc_tethered_fjc(L_two,p_m,0,s);
        V_two=log(kon(1)*c/koff(1));   % Unbinding energy of C2CD tethered by C2E in the absence of force

        Vub=log(kmax./koff);   % Unbinding energy barrier in kT.

        % Energy landscape in the absence of force
        Vm=[-V_two-V_one, -V_two-V_one+Vub(1), -V_one, -V_one+Vub(2) 0];

        r=(D-Ls)./Lm;  % Extension to contour length ratio
        r(r<0)=0;
        r(end)=0;  % The unbound state is not stretched.
        [force, em]=MS(r,p_m);   % Stretching force and energy per contour length
        em=Lm.*em./4.1;  % Entropic energy of the unfolded polypeptide
        energy=em+Vm+V;
        e=energy(1:2:end);   % Stable state energy
        fp=force(1:2:end);  % State pulling force
        
end

