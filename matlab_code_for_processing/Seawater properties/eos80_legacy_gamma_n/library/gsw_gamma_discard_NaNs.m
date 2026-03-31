function [salinity_bottle,temperature_bottle,pressure_bottle,salinity_cast,temperature_cast,pressure_cast,Idata] = gsw_gamma_discard_NaNs(salinity_bottle,temperature_bottle,pressure_bottle,salinity_cast,temperature_cast,pressure_cast,Idata)

% gsw_gamma_discard_NaNs                  eliminate NaN's from the data set
%==========================================================================
%
% USAGE:
%  [salinity_bottle, temperature_bottle, pressure_bottle,...
%   salinity_cast, temperature_cast, pressure_cast, Idata] = ...
%   gsw_gamma_discard_NaNs(salinity_bottle,temperature_bottle,pressure_bottle,...
%                       salinity_cast,temperature_cast,pressure_cast,Idata)
%
% DESCRIPTION:
%  This programme removes NaN's from the data. The function is a subroutine
%  of the neatural density code.
%
% INPUT:
%  data     =  data values
%
% Optional
%  weights  =  weights of the data
%
% OUTPUT:
%  weighted_mean     =  weighted mean of the data
%  weighted_weights  =  weighted weights of the weights
%
% AUTHOR: David Jackett
%  Modified by Guillaume Serazin, Stefan Riha and Paul Barker
%
% VERSION NUMBER: 3.05.9 (28th January, 2017)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Jackett, D.R., and T.J. McDougall, 1997: A Neutral Density Variable for 
%   the World's Oceans. J. Phys. Oceanogr., 27, pp. 237–263.
%
%  McDougall, T.J., 1987: The vertical motion of submesoscale coherent
%   vortices across neutral surfaces. Journal of physical oceanography, 
%   17, pp.2334-2342.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

Iwet = ~isnan(salinity_bottle);
Idry = all(isnan(salinity_cast));
Iwet = Iwet & ~Idry;

salinity_bottle = salinity_bottle(Iwet);
temperature_bottle = temperature_bottle(Iwet);
pressure_bottle = pressure_bottle(Iwet);

salinity_cast = salinity_cast(:,Iwet);
temperature_cast = temperature_cast(:,Iwet);
pressure_cast = pressure_cast(:,Iwet);

Idata = Idata(Iwet);

end