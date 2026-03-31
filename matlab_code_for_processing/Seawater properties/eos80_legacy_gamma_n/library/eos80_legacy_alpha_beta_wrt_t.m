function [alpha_wrt_t, beta_const_t] = eos80_legacy_alpha_beta_wrt_t(SP,t,p)

% eos80_legacy_alpha_beta_wrt_t               thermal expansion coefficient 
%                            with respect to in-situ temperature and saline 
%                   contraction coefficient at constant in-situ temperature
%                                                                   (eos80)
%==========================================================================
%
% USAGE:  
%  alpha_wrt_t, beta_const_t = eos80_legacy_alpha_beta_wrt_t(SA,t,p)
%
% DESCRIPTION:
%  Calculates the thermal expansion coefficient of seawater with respect to  
%  in-situ temperature and the saline (i.e. haline) contraction coefficient
%  of seawater at constant in-situ temperature.
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
%  alpha_wrt_t  =  thermal expansion coefficient                    [ 1/K ]
%                  with respect to in-situ temperature
%  beta_const_t =  saline contraction coefficient              [ unitless ]
%                  at constant in-situ temperature
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
%  Millero, et. al., 1980: Equation of state for seawater. Deep-Sea Res, 
%   27a, pp 255-264.
%
%  Millero and Poisson, 1981: Deep-Sea Res.,28A PP 625-629.
%   both of the baove references are also found in UNESCO report 38 (1981)
%
%==========================================================================

p_scaled = p.*0.1;
sqrt_SP = sqrt(SP);

r3500 = 1028.1063;
r4 = 4.8314e-4;
dr350 = 28.106331;

% *********************************************************
% pure water density at atmospheric pressure
% Bigg, p.h., 1967: br. J. Applied Physics, 8, pp 521-537.
% *********************************************************
r1 = ((((6.536332e-9.*t - 1.120083e-6).*t + 1.001685e-4).*t - 9.095290e-3).*t + 6.793952e-2).*t - 28.263737;

% seawater density atm press.
% coefficients involving salinity
r2 = (((5.3875e-9.*t - 8.2467e-7).*t + 7.6438e-5).*t - 4.0899e-3).*t + 8.24493e-1;
r3 = (-1.6546e-6.*t + 1.0227e-4).*t - 5.72466e-3;

%  international one-atmosphere equation of state of seawater
sig = (r4.*SP + r3.*sqrt_SP + r2).*SP + r1;

% specific volume at atmospheric pressure
v350p0 = 9.726620681149411e-04; %v350p0 = 1./r3500;
sva_p0 = -sig.*v350p0./(r3500 + sig);

sigma_p0 = sig + dr350;

spec_vol_p0 = 1./(1000 + sigma_p0);

% Calculate the derivative of rho wrt SP
r4s = 9.6628e-4;
rho_SP = r4s.*SP + 1.5.*r3.*sqrt_SP + r2;

% Calculate the derivative of rho wrt t
r1 =(((3.268166e-8.*t - 4.480332e-6).*t + 3.005055e-4).*t - 1.819058e-2).*t + 6.793952e-2;
r2 = ((2.155e-8.*t - 2.47401e-6).*t + 1.52876e-4).*t - 4.0899e-3;
r3 = -3.3092e-6.*t + 1.0227e-4;
rho_p0_t = (r3.*sqrt_SP + r2).*SP + r1;

rho_p0 = 1000 + sigma_p0;
rho_p02 = rho_p0.*rho_p0;

specvol_p0_t = -rho_p0_t./rho_p02;

% Calculate the derivative of specific volume wrt SP
specvol_SP = -rho_SP./rho_p02;

% ******************************************************************
% Equation of state for seawater
% Millero, et. al., 1980: DSR 27a, pp 255-264.
% ********************************************************
% compute the compression terms
b_SP = (9.1697e-10.*t + 2.0816e-8).*t - 9.9348e-7;
bw = (5.2787e-8.*t - 6.12293e-6).*t + 3.47718e-5;
b = bw + b_SP.*SP;

% Calculate the derivative of beta wrt t
bw = 1.05574e-7.*t - 6.12293e-6;
e = 1.83394e-9.*t + 2.0816e-8;
b_t = bw + e.*SP;

d = 1.91075e-4;
c = (-1.6078e-6.*t - 1.0981e-5).*t + 2.2838e-3;
aw = ((-5.77905e-7.*t + 1.16092e-4).*t + 1.43713e-3).*t - 0.1194975;
a = (d.*sqrt_SP + c).*SP + aw;

% Calculate the derivative of a wrt SP
a_SP = 2.866125e-4.*sqrt_SP + c;

% Calculate the derivative of a
c = -3.2156e-6.*t - 1.0981e-5;
aw = (-1.733715e-6.*t + 2.32184e-4).*t + 1.43713e-3;
a_t = c.*SP + aw;

% Calculate the COEFFICIENT k0
b1 = (-5.3009e-4.*t + 1.6483e-2).*t + 7.944e-2;
a1 = ((-6.1670e-5.*t + 1.09987e-2).*t - 0.603459).*t + 54.6746;
kw = (((-5.155288e-5.*t + 1.360477e-2).*t - 2.327105).*t + 148.4206).*t - 1930.06;
k0 = (b1.*sqrt_SP + a1).*SP + kw;

% Calculate the derivative of k0 wrt SP
k0_SP = 1.5.*b1.*sqrt_SP + a1;

% Calculate the derivative of k wrt SP
k_SP = (b_SP.*p_scaled + a_SP).*p_scaled + k0_SP;

% Calculate the derivative of k0 wrt t
b1 = -1.06018e-3.*t + 1.6483e-2;
a1 = (-1.8501e-4.*t + 2.19974e-2).*t - 0.603459;
kw = ((-2.0621152e-4.*t + 4.081431e-2).*t - 4.65421).*t + 148.4206;
k0_t = (b1.*sqrt_SP + a1).*SP + kw;

% EVALUATE PRESSURE POLYNOMIAL
% ***********************************************
%   K EQUALS THE SECANT BULK MODULUS OF SEAWATER
%   DK = K(S,T,P)-K(35,0,P)
%   K35 = K(35,0,P)
% ***********************************************
dk = (b.*p_scaled + a).*p_scaled + k0;
k35  = (5.03217e-5.*p_scaled + 3.359406).*p_scaled + 21582.27;
gam = p_scaled./k35;
pk = 1 - gam;
sva = sva_p0.*pk + (v350p0 + sva_p0).*p_scaled.*dk./(k35.*(k35 + dk));

v350p = v350p0.*pk;

%  ****************************************************
% COMPUTE DENSITY ANAMOLY WITH RESPECT TO 1000.0 KG/M**3
%  1) DR350: DENSITY ANAMOLY AT 35 (IPSS-78), 0 DEG. C AND 0 DECIBARS
%  2) DR35P: DENSITY ANAMOLY 35 (IPSS-78), 0 DEG. C ,  PRES. VARIATION
%  3) DVAN : DENSITY ANAMOLY VARIATIONS INVOLVING SPECFIC VOL. ANAMOLY
% ********************************************************************
% CHECK VALUE: SIGMA = 59.82037  KG/M**3 FOR S = 40 (IPSS-78),
% T = 40 DEG C, P0= 10000 DECIBARS.
% *******************************************************

dr35p = gam./v350p;
dvan = sva./(v350p.*(v350p + sva));
sigma = dr350 + dr35p - dvan;

k = k35 + dk;
rec_k = 1./k;
rec_k2 = rec_k.*rec_k;

specvol_p = 1 - p_scaled.*rec_k;
k_t = (b_t.*p_scaled + a_t).*p_scaled + k0_t;

rho = sigma + 1000;
specvol = 1./rho;
rho2 = rho.*rho; % 1./(specvol.*specvol);

% Calculate the derivative of specvol wrt SP
specvol_SP = specvol_SP.*specvol_p + spec_vol_p0.*p_scaled.*k_SP.*rec_k2;
rho_SP = -specvol_SP.*rho2;
beta_const_t = rho_SP.*specvol;

specvol_t = specvol_p0_t.*specvol_p + spec_vol_p0.*p_scaled.*k_t.*rec_k2;
rho_t = -specvol_t.*rho2;
alpha_wrt_t = rho_t.*specvol;

end


