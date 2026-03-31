function [o2s,KHo2]=gassat07_o2only(temp,sal)
% [n2s,o2s,ars,nes,KHn2,KHar,KHne,KHo2]=gassat07_o2only(temp,sal)

% GASSAT (Jan 2, 2004) Corrected Ne, Ar and N2 on June 2007
% function program to calculate saturation concentrations and Henry's law
% coefficients for N2,O2,Ar,CO2,Ne,He and H2O as a function of temp (C)
% and sal (ppt)
% Calculations are done at 1 atmosphere in umol/kg from the best available sources
% in 2004:            
%	O2; Garcia and Gordon, 1992, L&O, 37, 1307 (reinterpretation of Benson and
%      Krause, 1984, L&O,29,662)
%  Ar, N2, and Ne; Hamme and Emerson, 2004, DSR I, 51, 1517
%  CO2; Weiss,1974, Mar. Chem. 2, 203
%  He; Weiss, 1971, J.Chem.Eng.Data,16,235 
%  H2O  DOE,1974,Dickson and Goyet,ORNL/CDIAC-74
%  Sea water density Millero and Poisson,1981,Deep-Sea Res.,28,625 as reported
%       in DOE,1974

tabs = temp + 273.15;
P = 1.0;
% Xn2 = .78084;
Xo2 = .20946;
% Xar = .00934;
% Xne = .0000182;
% Xhe = 5.24E-06;
% Xco2 = 3.65E-04;

% MOLAR VOLUMES (liters / mole; don't need it for CO2 because KH,CO2 is in molar units)

% molvolo2 = 22.3916;
% molvoln2 = 22.3919;
% molvolar = 22.3869;
% molvolne = 22.4224;
% molvolhe = 22.1414;

% CALCULATE H&E N2s,ARs,Nes; G&G O2; Weiss CO2s and He in mol/kg 
% for CO2 CALCULATE, CO2s,FROM THE HENRY'S LAW COEFFICIENT

o2s = GGo2(temp,sal)*1e-06; %mol/kg
% ars = HEar(temp,sal)*1e-06;
% n2s = HEn2(temp,sal)*1e-06;
% nes = HEne(temp,sal)*1e-09;
% co2s =weissKHco2(tabs,sal)*Xco2;
% densw = sigma(temp,sal);
% hes = weissHe(tabs,sal)/(molvolhe*densw);
% n2ar = n2s/ars; 

% CALUCALTE HENRY'S LAW COEFFICIENTS, KH (mol kg-1 atm-1), FROM SATURATION CONCENTRATION, MOLE FRACTION,
% Xi, AND THE SATURATION PRESSURE OF WATER, ph2o.

ph2os = h2osat(tabs,sal);
delp = P-ph2os;
KHo2 = o2s ./ (Xo2.*delp);
% KHn2 = n2s / (Xn2*delp);
% KHar = ars / (Xar*delp);
% KHne = nes / (Xne*delp);
% KHhe = hes / (Xhe*delp);
% KHco2=weissKHco2(tabs,sal);

function [o2s]=GGo2(T, S)

% GGO2 This function calculates the oxygen saturation concentration in umol/kg using the 
%  regression of Garcia and Gordon, 1992,L&O,1307 with Benson and Krause
%  data

TS = log((298.15 - T) ./ (273.15 + T));
A0 = 5.80871; A1 = 3.20291; A2 = 4.17887; A3 = 5.10006; A4 = -.0986643; A5 = 3.80369;
B0 = -.00701577; B1 = -.00770028; B2 = -.0113864; B3 = -.00951519; C0 = -2.75915E-07;
C = A0 + A1 .* TS + A2 .* TS .^ 2 + A3 .* TS .^ 3 + A4 .* TS .^ 4 + A5 .* TS .^ 5 + S .*...
   (B0 + B1 .* TS + B2 .* TS .^ 2 + B3 .* TS .^ 3) + C0 .* S .^ 2; 
o2s = exp(C);
return

function [ars]=HEar(T, S)
% HEar This function calculates the argon satruation value as a function of T and S
%         in mol/kg using the regression equations of Hamme and Emerson,
%         Mar. Chem., 2004, xx,xx

TS = log((298.15 - T) / (273.15 + T));
A0 = 2.7915; A1 = 3.17609; A2 = 4.13116; A3 = 4.90379; A4 = 0.0; A5 = 0.0;
B0 = -.00696233; B1 = -.00766673; B2 = -.0116888; B3 = 0.0; C0 = 0.0;
C = A0 + A1 * TS + A2 * TS ^ 2 + A3 * TS ^ 3 + A4 * TS ^ 4 + A5 * TS ^ 5 + S *...
   (B0 + B1 * TS + B2 * TS ^ 2 + B3 * TS ^ 3) + C0 * S ^ 2; 
ars = exp(C);
return

function [n2s]=HEn2 (T, S)
% HEn2 This function calculates the nitrogen satruation value as a function of T and S
%         in mol/kg using the regression equations of Hamme and Emerson,
%         Mar. Chem., 2004, xx,xx

TS = log((298.15 - T) / (273.15 + T));
A0 = 6.42931; A1 = 2.92704; A2 = 4.32531; A3 = 4.69149; A4 = 0.0; A5 = 0.0;
B0 = -.00744129; B1 = -.00802566; B2 = -.0146775; B3 = 0.0; C0 = 0.0;
C = A0 + A1 * TS + A2 * TS ^ 2 + A3 * TS ^ 3 + A4 * TS ^ 4 + A5 * TS ^ 5 + S *...
   (B0 + B1 * TS + B2 * TS ^ 2 + B3 * TS ^ 3) + C0 * S ^ 2; 
n2s = exp(C);
return

function [nes]=HEne(T,S)
% HEne This function calculates the neon satruation value as a function of T and S
%         in mol/kg using the regression equations of Hamme and Emerson,
%         Mar. Chem., 2004, xx,xx

TS = log((298.15 - T) / (273.15 + T));
A0 = 2.18156; A1 = 1.29108; A2 = 2.12504; A3 = 0.0; A4 = 0.0; A5 = 0.0;
B0 = -.00594737; B1 = -.00513896; B2 = 0.0; B3 = 0.0; C0 = 0.0;
C = A0 + A1 * TS + A2 * TS ^ 2 + A3 * TS ^ 3 + A4 * TS ^ 4 + A5 * TS ^ 5 + S *...
   (B0 + B1 * TS + B2 * TS ^ 2 + B3 * TS ^ 3) + C0 * S ^ 2; 
nes = exp(C);
return

function [KHco2] = weissKHco2(TT, S)
% WEISSCO2 This function calucaltes the Henry's Law coefficient for CO2 as a function
%          of T and S in mol/kg atm using the equations of Weiss,1974,Mar.Chem,2,203.

A1 = -60.2409; A2 = 93.4517; A3 = 23.3585;
B1 = .023517; B2 = -.023656; B3 = .0047036;
C = A1 + A2 * (100 / TT) + A3 * log(TT / 100) + S * (B1 + B2 * (TT / 100) + B3 * (TT / 100) ^ 2);
KHco2 = exp(C);
return

function [hes_ml]=weissHe(TT,S)
%	WEISSHe This function calculates the helium saturation value as a function of T and S
%          in m/l using the regression equations of Weiss,1971,J.Chem.Eng.,16,235. 
A1 = -152.9405; A2 = 196.8840; A3 = 126.8015; A4 = -20.6767;
B1 = -.040543; B2 = .021315; B3 = -.0030732;
C = A1 + A2 * (100 / TT) + A3 * log(TT / 100) + A4 * (TT / 100) + S * (B1 + B2 * (TT / 100) + B3 * (TT / 100) ^ 2);
hes_ml = exp(C);
return


function [ph2os]=h2osat(TT,S)
% H2OSAT This function calculates the saturation water vapor pressure in atm from 
%      Zeebe and Wolf-Gladrow (2000). There is a more accurate one in DOE (1994)

A1=24.4543; A2=6745.09; A3=4.8489; A4=0.000544;
C = A1 - A2./TT - A3*log(TT./100) - S.*A4;
ph2os = exp(C);
return

function densw = sigma(T,S)
% SIGMA This function calculates the density of seawater in kg/m3 as a function of 
%		  t and salinity. DOE,1994,Dickson and Goyet, ORNL/CDIAC-74

A1 = 0.824493; A2 = -.0040899; A3 = 7.6438E-05; A4 = -8.2467E-07; A5 = 5.3875E-09;
B1 = -.00572466; B2 = 1.0227E-04; B3 = -1.6546E-06;
C1 = 4.8314E-04;
D1 = 999.842594; D2 = .06793952; D3 = -0.009095290; D4 = .0001001685; D5 = -1.120083E-06;
D6 = 6.536332E-09;
SMOW = D1 + D2*T + D3*T^2 + D4*T^3 + D5*T^4 + D6*T^5;
A = A1 + A2*T + A3*T^2 + A4*T^3 + A5*T^4;
B = B1 + B2*T + B3*T^2;
densw = SMOW + A*S + B*S^1.5 + C1*S^2;
return