function [ph2os]=h2osat_smb(T_C,S)
% H2OSAT This function calculates the saturation water vapor pressure in atm from 
%      Zeebe and Wolf-Gladrow (2000). There is a more accurate one in DOE (1994)
% copied from R. Hamme's gassat07.m and adapted to allow matrices.

TT = T_C + 273.15;
A1=24.4543; A2=6745.09; A3=4.8489; A4=0.000544;
C = A1 - A2./TT - A3.*log(TT./100) - S.*A4;
ph2os = exp(C);


end