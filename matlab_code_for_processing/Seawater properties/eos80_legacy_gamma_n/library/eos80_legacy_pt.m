function pt = eos80_legacy_pt(SP,t,p,pr)

% eos80_legacy_pt                             potential temperature (eos80)
% =========================================================================
%
% USAGE:
%  pt = eos80_legacy_pt(SP,t,p,p_ref)
%
% DESCRIPTION:
%  Calculates potential temperature with the general reference pressure, 
%  p_ref, from in-situ temperature, t, using Bryden (1973) polynomial for 
%  adiabatic lapse rate and runge-kutta 4-th order integration algorithm. 
%
% INPUT:
%  SP  =  Practical Salinity                                   [ unitless ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  p_ref  =  reference pressure                                    [ dbar ]
%
% OUTPUT:
%  pt  =  potential temperature with reference pressure, p_ref, on the 
%         ITS-90 temperature scale                                [ deg C ]
%
% AUTHOR: 
%  David Jackett, Paul Barker and Trevor McDougall     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05.9 (1st March, 2017)
%
%  checkvalue: 
%  SP = 40 (ipss-78), t = 40 deg c, p = 10000 decibars, pr = 0 decibars
%  pt = 36.89073 c.
%
% REFERENCES:            
%  Bryden, H., 1973: deep-sea res.,20,401-408.
%
%  Fofonoff, N., 1977: deep-sea res.,24,489-491
%
%==========================================================================

[ms,ns] = size(SP);
[mp,np] = size(p);

if ms ~= mp & ms == np
    p = p.';
end

[mpr,npr] = size(pr);
if ms ~= mpr & ms == npr & ns == mpr 
    pr = pr.';
end

h = pr - p;
xk = h.*eos80_legacy_adiabatic_lapse_rate_from_t(SP,t,p);
t = t + 0.5.*xk;

q = xk;
p = p + 0.5.*h;
xk = h.*eos80_legacy_adiabatic_lapse_rate_from_t(SP,t,p);
t = t + 0.29289322.*(xk - q);

q = 0.58578644.*xk + 0.121320344.*q;
xk = h.*eos80_legacy_adiabatic_lapse_rate_from_t(SP,t,p);
t = t + 1.707106781.*(xk - q);

q = 3.414213562.*xk - 4.121320344.*q;
p = p + 0.5.*h;
xk = h.*eos80_legacy_adiabatic_lapse_rate_from_t(SP,t,p);
pt = t + (xk - 2.*q)./6;

end     