function [dn2,do2,dar,dco2,dne,Sn2,Sar]=gasdif(temp,sal)

% GASDIF (feb, 2001)
% function program to calculate the molecular diffucion coefficents,d, and
% and schmidt numbers, S, of gases as a funciton of temperature, t (deg C)
% for fresh water,sal=0, and seawater, sal = 35 (ppt)
% The Eyring (Arrhenius) equation is used for diffucion constants
%     preexponential factors, A, and acitivation energies, E, are from:
%		for N2, O2, and Ar:  Wise and Houghton, 1966, Chem Eng. Sci,21, 994
%		for the rest; Jahne et al., 1987, JGR, 92, 10,767
% The kinematic viscosity for fresh water and seawater are fit to the data
%     of Sverdrup,Johnson and Flemming (The Oceans, 1942,Wiley Interscience)
%		between 5 and 25 deg C (not including the data at 0 and 35 deg)using
%		an empirical exponential function

tabs = temp + 273.15;
R = 0.008314;

% Pre-exponential terms (cm2 s-1)

An2=0.07900; Ao2=.04200; Aar=0.10600; Aco2=0.05019; Ane=0.01608; Ahe=0.00818;
Akr=0.06393; Axe=0.09007; Ach4=0.03047;

% Activation Energies (kJoul mol-1)

En2=19.62; Eo2=18.37; Ear=20.63; Eco2=19.51; Ene=14.84; Ehe= 11.70; Ekr=20.20;
Exe=21.61; Ech4=18.36;

% Kinematic viscosity at sal = 0

lnv= -.02662*temp - 4.062;
v=exp(lnv);

% Changes if sal = 35

C=1;
if sal==35.0
   C=0.94;
   lnv = -.02580*temp - 4.036;
   v = exp(lnv);
end

% Calculate molecular diffusion coefficients (dx, cm2 s-1)

RT=R*tabs;
dn2=C*An2*exp(-En2/RT); do2=C*Ao2*exp(-Eo2/RT); dar=C*Aar*exp(-Ear/RT);
dco2=C*Aco2*exp(-Eco2/RT); dne=C*Ane*exp(-Ene/RT); dhe=C*Ahe*exp(-Ehe/RT);
dkr=C*Akr*exp(-Ekr/RT); dxe=C*Axe*exp(-Exe/RT);dch4=C*Ach4*exp(-Ech4/RT);

% Calculate Schmid Numbers (S)

Sn2=v/dn2; So2=v/do2; Sar=v/dar; Sco2=v/dco2; Sne=v/dne; She=v/dhe; Skr=v/dkr;
Sxe=v/dxe; Sch4=v/dch4;

return