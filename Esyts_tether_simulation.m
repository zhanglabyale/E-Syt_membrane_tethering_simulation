% This script simulates trans-membrane binding of single monomers of E-Syt1 and E-Syt2. 
kT=4.1;
p_m=0.6;  % Persistence length of disordered polypeptides in nm

% Parameters for the membrane repulsive potential. See Eq. 6.
Emax=35; % Maximum membrane repulsive potential
dc=1;  
d1=10;

disp('Emax, dc, d1:')
disp([Emax, dc, d1])


leg={'E-Syt1','E-Syt2'};
kmax=1e6;  % Diffusion-limited binding rate
s=0.7;   % nm^2: surface area per lipid
kon=[4.8 4.6; 5.8 7.2];   % log10(kon), [C2CD C2E; C2AB C2C]
koff=[0.3 1.8; 2.6 1.5];  % log10(koff)
kon=10.^kon;     % kon (M-1s-1)
koff=10.^koff;   % koff (s-1)

Lmax=[191 181]; % Contour length of completely unfolded state in a.a. [E-Syt1 E-Syt2]
Ls=[4,6,6,8,8; 4,6,6,8,8];  % nm. Hard core sizes of C2 domains in the pulling direction [E-Syt1; E-Syt2]
Lm=[121 121 Lmax(1) Lmax(1) Lmax(1); 54 54 Lmax(2) Lmax(2) Lmax(2)];  % Contour length of the polypeptide bridging two membranes in a.a.
L_two=[70 127];   % Contour length of the tether for C2CD when C2E is bound to the trans-membrane.

% Contour length in nm
Lmax=Lmax*0.365;
Lm=Lm.*0.365;
L_two=L_two.*0.365;  

[c,e_cor] = conc_tethered_fjc(L_two,p_m,0,s);
V_two=log(kon(:,1)'.*c./koff(:,1)');   % Unbinding energy of C2CD tethered by C2E in the absence of force

Vub=log(kmax./koff);   % Unbinding energy barrier in kT.
Dmax=Lmax+Ls(3);  % Maximum tether length
N=100;  % Number of membrane distanc
D=linspace(5,35,N); 
M=length(Lm);  % Number of states and barriers

figure 
for j=1:2  % Different E-Syt: 1 for E-Syt1 and 2 for E-Syt2
    disp(leg{j})
    
    % Calculate equilibrium membrane separation,E-Syt tethering force, and
    % total energy of states with two or one C2 domain bound
    dmin=zeros(2,1);
    fm=zeros(2,1);
    ee=zeros(2,1);
    for k=1:2   % The two bound states: 1 for both C2 domains bound, 2 for only C-terminal C2 domain bound 
        m=2*k-1;
        % Find the equilibrium membrane distance where E-Syt tethering
        % force equates the membrane repulsive force
        fun = @(x)pullforce1(x,Lm(j,m),Ls(j,m),p_m,Emax,dc,d1);
        [dmin(k),value]=fzero(fun,15);
        if(abs(value)>1e-6)
           disp('Root not found')
           disp(value)
        end
        [~,fm(k),a] = force_energy2(dmin(k),kmax,kon(j,:),koff(j,:),Lm(j,:),Ls(j,:),p_m,s,Emax,dc,d1);
%         disp('State number, memebrane separation (nm), minumum state energy (kT)')
        ee(k)=a(k);
        disp([k dmin(k) ee(k) fm(k)])
    end
 
% force and energy of both stable and transition states for N membrane
% distances
    force=zeros(N,M);
    energy=zeros(N,M);
% Population of state states
    pop=zeros(N,3);
    pop_noCa=zeros(N,2);  % Population of C2E bound and unbound states in the absence of Ca
    lifetime=zeros(N,3);
    fav=zeros(N,1);   % Average force over all states in the presence of Ca
    rate=zeros(N,4);  % Unbinding and binding rate of C2CD and C2E
    fav_noCa=zeros(N,1); % Average force over all states in the absence of Ca

%   % Parameters for E-Syt1 without C2E
%    force1=zeros(N,2);
    energy1=zeros(N,2);
    pop1=zeros(N,2);
%     lifetime1=zeros(N,2);
     fav1=zeros(N,1);   % Average force
%     rate1=zeros(N,2);  % Unbinding and binding rate of C2CD and C2E

    f=zeros(N,1); % Membrane repulsive force

    for i=1:N   
        
%         [ff,ee] = force_energy2(D(i),kmax,kon(j,:),koff(j,:),Lm(j,:),Ls(j,:),p_m,s,Emax,dc,d1);
        
        % Membrane repusive potential and force 
        [V,f(i)]=membrane_potential_f2(D(i),Emax,dc,d1);
%         V=0;
%         f(i)=0;
        [c,e_cor] = conc_tethered_fjc(Lmax(j),p_m,0,s);
        kb=kon(j,2).*c;  
        V_one=log(kb/koff(j,2));  % Unbinding energy of C2E in kT at minimum membrane separation

        % Bound C2CD without C2E bound
        [c1,~] = conc_tethered_fjc(Lm(j,1),p_m,0,s);
        kb1=kon(j,1).*c1;  
        V_one1=log(kb1/koff(j,1));  % Unbinding energy of C2CD in kT at minimum membrane separation
        
        % Energy landscape in the absence of force
        Vm=[-V_two(j)-V_one, -V_two(j)-V_one+Vub(j,1), -V_one, -V_one+Vub(j,2) 0];

        r=(D(i)-Ls(j,:))./Lm(j,:);  % Extension to contour length ratio
        r(r<0)=0;
        r(end)=0;  % The unbound state is not stretched.
        [force(i,:), em]=MS(r,p_m);   % Stretching force and energy per contour length
%         em(em==Inf)=100;
        em=Lm(j,:).*em./kT;  % Entropic energy of the unfolded polypeptide
        energy(i,:)=em+Vm+V;
        e=energy(i,1:2:5);   % Stable state energy
        a=exp(-e);
        pop(i,:)=a./sum(a);  % State population
        a_noCa=[a(2) a(3)];
        pop_noCa(i,:)=a_noCa./sum(a_noCa);   
        fav(i)=sum(force(i,1:2:5).*pop(i,:));  % Average net force
        fav_noCa(i)=sum(force(i,3:2:5).*pop_noCa(i,:));
       
        energy1(i,:)=[em(1)-V_one1 0]+V;   % C2CD binding in the absence of C2E
        a=exp(-energy1(i,:));
        pop1(i,:)=a./sum(a);
        fav1(i)=force(i,1).*pop1(i,1);
        
        rate(i,1:2:3)=kmax.*exp(energy(i,1:2:3)-energy(i,2:2:4));  % Unbinding rate
        rate(i,2:2:4)=kmax.*exp(energy(i,3:2:5)-energy(i,2:2:4));  % Binding rate
    %     de(2,2)=log(kmax/kb);   % Binding energy barrier of C2E
        RATE=rate(i,:);
        RATE(RATE>1e8)=NaN;
        rate(i,:)=RATE;
        lifetime(i,:)=[1./rate(i,1) 1./(rate(i,2)+rate(i,3)) 1./rate(i,4)];

    end

    d=D;

    a=0.8/3;
    b=(j-1)*0.5;

    ax1=subplot('position',[0.1+b 0.1+2*a 0.35 a]);   
    plot(d,pop,d,pop_noCa,d,pop1)
    set(gca,'box','on','xticklabel','')
    grid on
    ylabel('Population')

    ax2=subplot('position',[0.1+b 0.1+a 0.35 a]); 
    % plot(d,force(:,1:2:end),d,fav,d,fav1)
    plot(d,force(:,1:2:end),d,f,d,fav,d,fav_noCa,d,fav1)
    set(gca,'box','on','xticklabel','')
    grid on
    ylabel('Force (pN)')
    ylim([0 10])
    hold on
     plot(dmin,fm,'ro')
    hold off

    ax3=subplot('position',[0.1+b 0.1 0.35 a]); 

    % plot(d,energy(:,1:2:end),d,energy1(:,1))
    plot(d,energy(:,1:2:end),d,energy1)
    set(gca,'box','on')
    grid on
    ylabel('Energy (kT)')
    ylim([-5 15])
    hold on
    plot(dmin,ee,'ro')
    hold off
    xlabel('Membrane distance (nm)')   
%     ax4=subplot('position',[0.1 0.1 0.8 a]); 
%     % semilogy(d,lifetime,d,lifetime1)
%     semilogy(d,lifetime)
%     set(gca,'box','on')
%     grid on
%     ylabel('Lifetime (s)')
%     xlabel('Membrane distance (nm)')

    linkaxes([ax1 ax2 ax3],'x')
    xlim([5 35])

end