function pot_rho_pref = eos80_legacy_pot_rho_pref(SP,t,p,pref)

% eos80_legacy_pot_rho_pref                               potential density
%==========================================================================
%
% USAGE:
%  pot_rho_pref = eos80_legacy_pot_rho_pref(SP,t,p,p_ref)
%
% DESCRIPTION:
%  Calculates potential density of seawater.  Note. This function outputs
%  potential density, not potential density anomaly; that is, 1000 kg/m^3
%  is not subtracted.  
%
% INPUT:
%  SP  =  Practical Salinity                                   [ unitless ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  p_ref  =  reference pressure                                    [ dbar ]
%            ( i.e. reference absolute pressure - 10.1325 dbar )
%
% OUTPUT:
%  pot_rho_pref  =  potential density (not potential density anomaly)
%                                                                [ kg/m^3 ]
% AUTHOR: 
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
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

pt = eos80_legacy_pt(SP,t,p,pref);

pot_rho_pref = 1000 + eos80_legacy_sigma(SP,pt,pref);

end     