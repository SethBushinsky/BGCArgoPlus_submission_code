function [ph2os]=ph2osat_smb(T_C,S)
% H2OSAT This function calculates the saturation water vapor pressure in atm from
%      Zeebe and Wolf-Gladrow (2001), p. 281. There is a more accurate one
%      in DOE (1994) - the difference between the Zeebe/DOE calculations
%      are ~0.1% (out of water vapor pressure, not total pressure).
% copied from R. Hamme's gassat07.m and adapted to allow matrices.
% Equation and coefficients are originally from Weiss and Price 1980.  In
% this reference, it also states that pH2O result is in atm, not mole
% fraction.  This is counter to Sarmiento and Gruber, p. 78, Panel 3.2.1.
% Total difference would be on the order of a few percent of water vapor
% pressure, so likely not critical.  
% Seth Bushinsky, smb4@hawaii.edu

% ZW = 1;

TT = T_C + 273.15;

% if ZW==1, %use Zeebe and Wolf-Gladrow(2001)
    
    
    A1=24.4543; A2=6745.09; A3=4.8489; A4=0.000544;
    C = A1 - A2./TT - A3.*log(TT./100) - S.*A4;
    ph2os = exp(C);
    
% elseif ZW==0,% use DOE 1994
%     
%     a0= 2794.0144;
%     a(1)= 1430.6181;
%     a(2)= -18.2465;
%     a(3)= 7.6875;
%     a(4)=-0.0328;
%     a(5)= 0.2728;
%     a(6)= 0.1371;
%     a(7)= 0.0629;
%     a(8)= 0.0261;
%     a(9)= 0.0200;
%     a(10)= 0.0117;
%     a(11)= 0.0067;
%     
%     x = (2.*(TT)-921)./375;
%     
%     sum_akEk = 0;
%     for k=1:11,
%         sum_akEk = a(k).*x + sum_akEk;
%     end
%     
%     ph2os = 10^((.5.*a0+sum_akEk)./TT);
%     
% end
end