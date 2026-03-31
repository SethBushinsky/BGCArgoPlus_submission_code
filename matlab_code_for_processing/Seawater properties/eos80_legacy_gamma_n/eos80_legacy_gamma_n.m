function [gamma_n, gamma_error_lower, gamma_error_upper] = eos80_legacy_gamma_n(SP,t,p,long,lat)

% eos80_legacy_gamma_n                                      neutral density
%==========================================================================
%
% USAGE:  
%  [gamma_n,{gamma_error_lower,gamma_error_upper}] = ...
%                             eos80_legacy_gamma_n(SP,t,p,long,lat)
%
% DESCRIPTION:
%  Calculates neutral density value, gamma_n, in the open ocean by 
%  spatially interpolating the global reference data set, gamma_n_ref, to
%  the location of the seawater sample.  The input variables are those of
%  EOS-80, namely Practical Salinity and in-situ temperature.
% 
%  This function uses the gamma_n look up table of Jackett & McDougall 
%  (1997). 
%
% INPUT:
%  SP    =  Practical Salinity                                 [ unitless ]
%  t     =  in-situ temperature (ITS-90)                          [ deg C ]
%  p     =  sea pressure                                           [ dbar ] 
%          ( i.e. absolute pressure - 10.1325 dbar )
%  long  =  longitude in decimal degrees                     [ 0 ... +360 ]
%                                                      or [ -180 ... +180 ]
%  lat   =  latitude in decimal degrees north               [ -90 ... +90 ]
%
%  SP and t need to have the same dimensions.
%  p, lat & long may have dimensions 1x1 or Mx1 or 1xN or MxN,
%  where SP & t is MxN.
%
% OUTPUT:
%  gamma_n  =  neutral density                                     [ g/kg ]
% Optional:
%  gamma_error_lower = lower estimate of the error
%  gamma_error_upper = upper estimate of the error
%
% AUTHOR: 
%  David Jackett                                       [ help@teos-10.org ]
%
% MODIFIED:
%  Paul Barker and Trevor McDougall
%
% VERSION NUMBER: 3.06 (19th September, 2017)
%
% REFERENCES:
% Jackett, D.R., and T.J. McDougall, 1997: A Neutral Density Variable for
%  the World's Oceans. J. Phys. Oceanogr., 27, 237–263. 
%  doi: 10.1175/1520-0485(1997)0272.0.CO;2
%
% Unesco, 1983: Algorithms for computation of fundamental properties of 
%  seawater. Unesco technical papers in marine science 44, 53pp.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 5)
   error('eos80_legacy_gamma_n:  Requires five inputs')
end

sz_SP = size(SP);
sz_t = size(t);

if sum(sz_SP - sz_t) ~= 0 
    error('eos80_legacy_gamma_n: SP and t need be of the same size')
end 

p = gsw_resize(p,sz_SP);
long = gsw_resize(long,sz_SP);
lat = gsw_resize(lat,sz_SP);

long(long < 0) = long(long < 0) + 360; 

if any(p < -1.5)
    error('eos80_legacy_gamma_n: pressure needs to be positive')
end

%set any pressures between 0 and 1.5 to be equal to 0 (i.e. the surface)
p(p < 0) = 0;

if any(long < 0 | long > 360 | lat < -80 | lat > 64)
    [Ioor] = find(long < 0 | long > 360 | lat < -80 | lat > 64);
    SP(Ioor) = NaN;
end

if any(SP(:) < 0 | SP(:) > 42 | t(:) < -2.5 | t(:) > 40 | p(:) < 0 | p(:) > 10000)
    [Ioor] = find(SP(:) < 0 | SP(:) > 42 | t(:) < -2.5 | t(:) > 40 | p(:) < 0 | p(:) > 10000);
    SP(Ioor) = NaN;
end

[Iocean] = find(~isnan(SP(:) + t(:) + p(:) + lat(:) + long(:)));

gamma_n = nan(size(SP));
if nargout > 1
    gamma_error_lower = gamma_n;
    gamma_error_upper = gamma_n;
    [gamma_n(Iocean), gamma_error_lower(Iocean), gamma_error_upper(Iocean)] = ...
        eos80_legacy_gamma_n_labelling(SP(Iocean),t(Iocean),p(Iocean),long(Iocean),lat(Iocean));
else
    gamma_n(Iocean) = eos80_legacy_gamma_n_labelling(SP(Iocean),t(Iocean),p(Iocean),long(Iocean),lat(Iocean));
end


end
