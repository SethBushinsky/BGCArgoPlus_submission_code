function [gamma, weight] = eos80_legacy_gamma_n_labelling_extrapolation(SP,t,p,SP_ref,t_ref,p_ref,gamma_n_ref)

% eos80_legacy_gamma_n_labelling_extrapolation                 extrapolates
%               shallower and deeper than the neutral density look-up table
%==========================================================================
%
% USAGE:
%  [gamma, weight] = eos80_legacy_gamma_n_labelling_extrapolation(SP,t,p,...
%                                     SP_ref,t_ref,p_ref,gamma_n_ref)
%
% DESCRIPTION:
%  This function extrapolates shallower and deeper than the neutral density
%  look-up table
%
%  The correct way of calculating gamma_n is by calling
%  eos80_legacy_gamma_n.
%
%  This function uses version 1.0 of the gamma_n look up table (1997).
%
% INPUT:
%  SP        = Practical Salinity of the bottle                [ unitless ]
%  t         = in-situ temperature (ITS-90) of the bottle         [ deg C ]
%  p         = sea pressure of the bottle                          [ dbar ]
%          ( i.e. absolute pressure - 10.1325 dbar )
%  SP_ref    = Practical Salinity of the reference cast        [ unitless ]
%  t_ref     = in-situ temperature of the reference cast          [ deg C ]
%  p_ref     = sea pressure of the reference cast                  [ dbar ]
%  gamma_n_ref = neutral density of the reference cast           [ kg/m^3 ]
%
%  SP, t & p to be vectors and have the same dimensions, 1xM.
%  SP_ref, t_ref, p_ref & gamma_n_ref to be vectors and have the same
%  dimensions. They need to have dimensions NxM. wher ethe first dimension
%  is the vertical dimension.
%
% OUTPUT:
%  gamma  =  neutral density                                       [ g/kg ]
%  weight =  weight, depending on the magnitude of the extrapolation
%
% AUTHOR:
%  David Jackett                                       [ help@teos-10.org ]
%
% MODIFIED:
%  Paul Barker and Trevor McDougall
%
% VERSION NUMBER: 3.05.9 (1st March, 2017)
%
% REFERENCES:
%  Jackett, D.R., and T.J. McDougall, 1997: A Neutral Density Variable for
%   the World's Oceans. J. Phys. Oceanogr., 27, pp. 237–263.
%
% Unesco, 1983: Algorithms for computation of fundamental properties of
%  seawater. Unesco technical papers in marine science 44, 53pp.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

gamma = nan(size(SP));
weight = ones(size(SP));

SP = SP(:).';
t = t(:).';
p = p(:).';

pmid_top = 0.5*(p + p_ref(1,:));

pot_rho_top = eos80_legacy_pot_rho_pref(SP,t,p,pmid_top);
pot_rho_ref_top = eos80_legacy_pot_rho_pref(SP_ref(1,:),t_ref(1,:),p_ref(1,:),pmid_top);

if any(pot_rho_top < pot_rho_ref_top)
    
    [Ilight] = find(pot_rho_top < pot_rho_ref_top);
    
    pmid_ref_top = 0.5*(p_ref(1) + p_ref(2));
    
    pot_rho_ref_shallow = eos80_legacy_pot_rho_pref(SP_ref(1,Ilight),t_ref(1,Ilight),p_ref(1),pmid_ref_top);
    pot_rho_ref_deep = eos80_legacy_pot_rho_pref(SP_ref(2,Ilight),t_ref(2,Ilight),p_ref(2),pmid_ref_top);
    
    b_ref = (gamma_n_ref(1,Ilight) - gamma_n_ref(2,Ilight))./(pot_rho_ref_shallow - pot_rho_ref_deep);
    
    gamma(Ilight) = gamma_n_ref(1,Ilight) + b_ref.*(pot_rho_top(Ilight) - pot_rho_ref_top(Ilight));
    
end

if any(isnan(gamma))
    
    [Inan] = find(isnan(gamma));
    Inan = Inan(:)';
    
    [dummy, Idata] = nanmax(gamma_n_ref(:,Inan));
    Idata = Idata(:)';
    
    pmid_bottom = 0.5*(p(Inan) + p_ref(Idata).');
    
    pot_rho_bottom = eos80_legacy_pot_rho_pref(SP(Inan),t(Inan),p(Inan),pmid_bottom);
    
    dummy_I = sub2ind(size(SP_ref),Idata,Inan);
    pot_rho_ref_bottom = eos80_legacy_pot_rho_pref(SP_ref(dummy_I),t_ref(dummy_I),p_ref(Idata).',pmid_bottom);
    
    if any(pot_rho_bottom >  pot_rho_ref_bottom)
        
        [Iheavy] = find(pot_rho_bottom >  pot_rho_ref_bottom);
        
        pmid_ref_bottom = 0.5*(p_ref(Idata(Iheavy)) + p_ref(Idata(Iheavy)-1));
        
        dummy_Ishallow = sub2ind(size(SP_ref),Idata(Iheavy)-1,Inan(Iheavy));
        
        pot_rho_ref_shallow = eos80_legacy_pot_rho_pref(SP_ref(dummy_Ishallow),t_ref(dummy_Ishallow),p_ref(Idata(Iheavy)-1),pmid_ref_bottom);
        pot_rho_ref_deep = eos80_legacy_pot_rho_pref(SP_ref(dummy_I(Iheavy)),t_ref(dummy_I(Iheavy)),p_ref(Idata(Iheavy)),pmid_ref_bottom);
        
        b_ref = (gamma_n_ref(dummy_I(Iheavy)) - gamma_n_ref(dummy_Ishallow))./(pot_rho_ref_deep - pot_rho_ref_shallow);
        
        delta_gamma =  b_ref.*(pot_rho_bottom(Iheavy) - pot_rho_ref_bottom(Iheavy));
        gamma(Inan(Iheavy)) = gamma_n_ref(dummy_I(Iheavy)) + delta_gamma;
        
        dummy = (0.3 - delta_gamma)./0.3;
        dummy(dummy < 0) = 0;
        weight(Inan(Iheavy)) = dummy;
        
    end
    
end

end