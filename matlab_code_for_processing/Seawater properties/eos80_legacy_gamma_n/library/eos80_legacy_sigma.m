function sigma = eos80_legacy_sigma(SP,t,p)

% eos80_legacy_sigma                             density of seawater - 1000 
%==========================================================================
%
% USAGE:  
%  sigma = eos80_legacy_sigma(SA,t,p)
%
% DESCRIPTION:
%  Calculates in-situ density - 1000 of seawater from Practical Salinity and 
%  in-situ temperature.  Note that the output, sigma, is density anomaly; 
%  that is, 1000 kg/m^3 is subracted from it.  
%
% INPUT:
%  SP  =  Practical Salinity                                   [ unitless ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  sigma  =  density anomaly                                     [ kg/m^3 ]
%    
% AUTHOR: 
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05.9 (1st March, 2017)
%
% REFERENCES:
%  Millero, et. al., 1980: Equation of state for seawater. DSR 27a, 
%   pp 255-264.
%  
%  Bigg, P.H., 1967: pure water density at atmospheric pressure.
%   J. Applied Physics, 8, pp 521-537.
%
%==========================================================================

SP = SP(:).';
t = t(:).';
p = p(:).';

p_scaled = p.*0.1;
sr = sqrt(SP);

r3500 = 1028.1063;
r4 = 4.8314e-4;
dr350 = 28.106331;

% *********************************************************
% pure water density at atmospheric pressure
% Bigg, p.h., 1967: br. J. Applied Physics, 8, pp 521-537.
% *********************************************************
r1 = ((((6.536332e-9.*t - 1.120083e-6).*t + 1.001685e-4).*t - 9.095290e-3).*t + 6.793952e-2).*t - 28.263737;

% seawater density atm press.
%  coefficients involving salinity
%  r2 = a   in notation of millero and poisson 1981
r2 = (((5.3875e-9.*t - 8.2467e-7).*t + 7.6438e-5).*t - 4.0899e-3).*t + 8.24493e-1;
%  r3 = b  in notation of millero and poisson 1981
r3 = (-1.6546e-6.*t + 1.0227e-4).*t - 5.72466e-3;
%  international one-atmosphere equation of state of seawater
sig = (r4.*SP + r3.*sr + r2).*SP + r1;

% specific volume at atmospheric pressure
v350p0 = 9.726620681149411e-04; %v350p0 = 1./r3500;
sva = -sig.*v350p0./(r3500 + sig);
sigma = sig + dr350;

if any(p_scaled > 0) 
    [Ip] = find(p_scaled > 0);
    % ******************************************************************
    % Equation of state for seawater 
    % Millero, et. al., 1980: DSR 27a, pp 255-264.
    % ********************************************************
    % compute compression terms
    e = (9.1697e-10.*t(Ip) + 2.0816e-8).*t(Ip) - 9.9348e-7;
    bw = (5.2787e-8.*t(Ip) - 6.12293e-6).*t(Ip) + 3.47718e-5;
    b = bw + e.*SP(Ip);
    
    d = 1.91075e-4;
    c = (-1.6078e-6.*t(Ip) - 1.0981e-5).*t(Ip) + 2.2838e-3;
    aw = ((-5.77905e-7*t(Ip) + 1.16092e-4).*t(Ip) + 1.43713e-3).*t(Ip) - 0.1194975;
    a = (d.*sr(Ip) + c).*SP(Ip) + aw;
    %
    b1 = (-5.3009e-4.*t(Ip) + 1.6483e-2).*t(Ip) + 7.944e-2;
    a1 = ((-6.1670e-5.*t(Ip) + 1.09987e-2).*t(Ip) - 0.603459).*t(Ip) + 54.6746;
    kw = (((-5.155288e-5.*t(Ip) + 1.360477e-2).*t(Ip) - 2.327105).*t(Ip) + 148.4206).*t(Ip) - 1930.06;
    k0 = (b1.*sr(Ip) + a1).*SP(Ip) + kw;
    
    % evaluate pressure polynomial
    % ***********************************************
    %   k equals the secant bulk modulus of seawater
    %   dk = k(SP,t,p) - k(35,0,p)
    %   k35 = k(35,0,p)
    % ***********************************************
    dk = (b.*p_scaled(Ip) + a).*p_scaled(Ip) + k0;
    k35  = (5.03217e-5.*p_scaled(Ip) + 3.359406).*p_scaled(Ip) + 21582.27;
    gam = p_scaled(Ip)./k35;
    pk = 1 - gam;
    sva(Ip) = sva(Ip).*pk + (v350p0 + sva(Ip)).*p_scaled(Ip).*dk./(k35.*(k35 + dk));
    
    v350p = v350p0.*pk;
    
    %  ****************************************************
    % compute density anamoly with respect to 1000.0 kg/m**3
    %  1) dr350: density anamoly at SP = 35 (ipss-78), t = 0 C and p = 0 decibars
    %  2) dr35p: density anamoly at SP = 35 (ipss-78), t = 0 C, and with pres. variation
    %  3) dvan : density anamoly variations involving specfic vol. anamoly
    % ********************************************************************
    % check value:
    % sigma(40,40,10000) = 59.82037 kg/m**3, (sigma(SP,t,p)) 
    % *******************************************************
    dr35p = gam./v350p;
    dvan = sva(Ip)./(v350p.*(v350p + sva(Ip)));
    sigma(Ip) = dr350 + dr35p - dvan;
end

end
