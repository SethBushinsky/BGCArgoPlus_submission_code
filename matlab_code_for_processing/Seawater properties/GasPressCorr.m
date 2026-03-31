function [correction_factor] = GasPressCorr(P,T,gas)

% GasPressCorr   Correction factor for hydrostatic pressure on gas partial pressure
%=========================================================================
% GasPressCorr Version 1.0; 7 Dec 2012
%          Author: Roberta C. Hamme (University of Victoria)
%
% USAGE:  correction_factor = GasPressCorr(P,T,gas)
%
% DESCRIPTION:
%    This correction factor is the ratio of the gas's partial pressure at 
%    the given hydrostatic pressure to that at the surface for an equal 
%    disssolved concentration
%
% INPUT:  (if P and T are not singular they must have same dimensions)
%   P = water pressure [dbar]
%   T = temperature [degree C]
%   gas = indicates which gas to make calculation for (1=N2, 2=O2, 3=Ar)
%
% OUTPUT:
%   correction_factor = ratio of partial pressures  [unitless] 
% 
% AUTHOR:  Roberta Hamme (rhamme@uvic.ca)
%
% REFERENCE:
%    equation based on Taylor, Craig D. (1978) The Effect of Pressure upon
%    the Solubility of Oxygen in Water, Archives of Biochemistry and Biophysics
%    191, pp. 375-384.
%
%    Partial Molal Volumes are my guess based on several papers, mainly
%    Enns et al. (1965) J. Phys. Chem., 69, 389-391
%    Taylor (1978) ref above
%    Moore et al. (1982) J. Chem. Eng. Data, 27, 22-24
%    Kennan and Pollack (1995) J. Chem. Phys., 93, 2724-2735
%    Zhou and Batino (2001) J. Chem. Eng. Data, 46, 331-332
%
% DISCLAIMER:
%    This software is provided "as is" without warranty of any kind.  
%=========================================================================

% CALLER: general purpose
% CALLEE: none

%----------------------
% Check input parameters
%----------------------
if nargin ~=3
   error('GasPressCorr.m: Must pass 3 parameters')
end %if

% CHECK S,T dimensions and verify consistent
[mp,np] = size(P);
[mt,nt] = size(T);
[mg,ng] = size(gas);

  
% Check that T&S have the same shape or are singular
if ((mp~=mt) || (np~=nt)) && (mp+np>2) && (mt+nt>2)
   error('GasPressCorr: P & T must have same dimensions or be singular')
end %if

% Check that gas is singular
if mg+ng > 2
   error('GasPressCorr: gas must be singular')
end

%------
% BEGIN
%------

% Define constants
R = 82.05746;          % Gas Constant in units mL atm mol^-1 K^-1

% Partial molal volume of the gas (mL mol^-1)
switch gas
    case 1  % N2
        Vc = 33;
    case 2  % O2
        Vc = 32;
    case 3  % Ar
        Vc = 32.5;
end

% Unit Conversions
P = P * 0.098692316931;     % convert pressure from dbar to atm
T = T + 273.15;             % convert temperature from celcius to Kelvin

% Make calculation
correction_factor = exp(Vc*P./(R*T));

