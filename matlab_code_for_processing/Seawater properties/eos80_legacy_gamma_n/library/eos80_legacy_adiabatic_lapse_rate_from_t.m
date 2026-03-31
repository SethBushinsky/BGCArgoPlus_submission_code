function adiabatic_lapse_rate = eos80_legacy_adiabatic_lapse_rate_from_t(SP,t,p) 

% eos80_legacy_adiabatic_lapse_rate          adiabatic temperature gradient
%                                                                   (eos80)
%==========================================================================
%
% USAGE:  
%  adiabatic = eos80_legacy_adiabatic_lapse_rate_from_t(SP,t,p)
%
% DESCRIPTION:
%  Calculates the adiabatic temperature gradient.
%
% INPUT:
%  SP  =  Practical Salinity                                   [ unitless ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  adiabatic_lapse_rate = adiabatic lapse rate               [ deg C/dbar ]
%
% AUTHOR: 
%  David Jackett, Paul Barker and Trevor McDougall     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05.9 (1st March, 2017)
%
% checkvalue: 
%  SP = 40 (ipss-78), t = 40 deg C, p = 10000 decibars
%  adiabatic_lapse_rate = 3.255976e-4 C/dbar
%
% REFERENCES:
%  Bryden, H., 1973: Deep-sea Res., 20, pp. 401-408
%
%==========================================================================

dSP = SP - 35;

adiabatic_lapse_rate = (((-2.1687e-16.*t + 1.8676e-14).*t - 4.6206e-13).*p ...
    + ((2.7759e-12.*t - 1.1351e-10).*dSP + ((-5.4481e-14.*t ...
    + 8.733e-12).*t - 6.7795e-10).*t + 1.8741e-8)).*p ...
    + (-4.2393e-8.*t + 1.8932e-6).*dSP ...
    + ((6.6228e-10.*t - 6.836e-8).*t + 8.5258e-6).*t + 3.5803e-5;


end