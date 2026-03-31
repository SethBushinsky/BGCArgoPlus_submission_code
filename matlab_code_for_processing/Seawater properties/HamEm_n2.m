function [n2s] = HamEm_n2 (T, S)
% outputs N2 in umnol kg-1
TS = log((298.15 - T) ./ (273.15 + T));


A0 = 6.42931 ;	
A1 = 2.92704 ;	
A2 = 4.32531 ;	
A3 = 4.69149 ;	
B0 = -0.00744129 ;	
B1 = -0.00802566 ;	
B2 = -0.0146775 ;

C = A0 + A1.*TS + A2.*TS.^2 + A3.*TS.^3 + S.*(B0 + B1.*TS + B2.*TS.^2);

n2s = exp(C);

end