function [alpha_wrt_pt, beta_const_pt] = eos80_legacy_alpha_beta_wrt_pt(SP,t,p)

% eos80_legacy_alpha_beta_wrt_pt               thermal expansion coefficient 
%                          with respect to potential temperature and saline 
%                 contraction coefficient at constant potential temperature
%                                                                   (eos80)
%==========================================================================
%
% USAGE:  
%  [alpha_wrt_pt, beta_const_pt] = eos80_legacy_alpha_beta_wrt_pt(SP,t,p)
%
% DESCRIPTION:
%  Calculates the thermal expansion coefficient of seawater with respect to  
%  potential temperature and the saline (i.e. haline) contraction
%  coefficient of seawater at constant potential temperature.
%
%  It is based on the following functions from Bob Millard 
%  - THETA(SP,t,p,p_ref);
%  - SVAN(SP,t,p,sigma);
%  - EOSED(SP,t,p,DRV)
%  while, alpha and beta are defined in terms of Practical Salinity and 
%  potential temperature (Gill, 1982, See SECTION 3.7.4).
%   
% INPUT:
%  SP  =  Practical Salinity                                   [ unitless ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  alpha_wrt_pt  =  thermal expansion coefficient                   [ 1/K ]
%                   with respect to potential temperature
%  beta_const_pt =  saline contraction coefficient             [ unitless ]
%                   at constant potential temperature
% AUTHOR: 
%  Trevor McDougall
%
% MODIFIED:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05.9 (1st March, 2017)
%
% REFERENCES: 
%  Bigg, P.H., 1967: pure water density at atmospheric pressure.
%   J. Applied Physics, 8, pp 521-537.
%
%  Gill, A.E., 1982: Atmosphere and Ocean Dynamics, Academic Press,
%   New York.
%
%  Millero, et. al., 1980: Equation of state for seawater. DSR 27a, 
%   pp 255-264.
%
%==========================================================================

[alpha_wrt_t, beta_const_t] = eos80_legacy_alpha_beta_wrt_t(SP,t,p);

A1 = -0.83198e-5;
A2 = 0.54065e-7;
A3 = -0.40274e-9;
B0 = -0.17439e-5;
B1 = 0.29778e-7;
C1 = 0.31628e-9;
C2 = -0.21987e-11;
D0 = 0.41057e-10;
E1 = -0.50484e-14;

dSP = SP - 35;

pt_t = 1 + p.*(A1 + t.*(2.*A2 + 3.*A3.*t) + dSP.*B1 + p.*(C1 + 2.*C2.*t + p.*E1));

pt_SP = p.*(B0 + B1.*t + p.*D0);

alpha_wrt_pt = alpha_wrt_t./pt_t;

beta_const_pt = beta_const_t + alpha_wrt_pt.*pt_SP;
     
end


