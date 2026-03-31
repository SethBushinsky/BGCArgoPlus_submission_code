function [SP_neutral_surface,t_neutral_surface,p_neutral_surface] = eos80_legacy_neutral_surfaces(SP,t,p,gamma_n,gamma_n_surfaces)

% eos80_legacy_neutral_surfaces                 Practical Salinity, in-situ
%                                   temperatures and pressures on specified
%                                         neutral density surfaces (EOS-80)
%==========================================================================
%
% USAGE:  
%  [SP_neutral_surface,t_neutral_surface,p_neutral_surface] = ...
%            eos80_legacy_neutral_surfaces(SP,t,p,gamma_n,gamma_n_surfaces)
%
% DESCRIPTION:
%  Calculates Practical Salinity, in-situ temperatures and pressures on
%  specified neutral density surfaces.  The input variables are those of
%  EOS-80, namely Practical Salinity and in-situ temperature.
%
% INPUT:
%  SP    =  Practical Salinity                                 [ unitless ]
%  t     =  in-situ temperature (ITS-90)                          [ deg C ]
%  p     =  sea pressure                                           [ dbar ] 
%          ( i.e. absolute pressure - 10.1325 dbar )
%  gamma_n  =  neutral density                                     [ g/kg ]
%  gamma_n_surfaces = specified neutral density surfaces
%
%  SP, t and gamma_n need to have the same dimensions.
%  p, lat & long may have dimensions 1x1 or Mx1 or 1xN or MxN,
%  where SP & t is MxN.
%  gamma_n_surfaces needs to be a vector.
%
% OUTPUT:
%  SP_neutral_surface = Practical Salinity on the neutral density surfaces
%  t_neutral_surface  = in-situ temperature on the neutral density surfaces
%  p_neutral_surface  = pressure on the neutral density surfaces
%
% AUTHOR: 
%  David Jackett                                       [ help@teos-10.org ]
%
% MODIFIED:
%  Paul Barker and Trevor McDougall
%
% VERSION NUMBER: 3.05.9 (9th March, 2017)
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
   error('eos80_legacy_neutral_surfaces:  Requires five inputs')
end

sz_SP = size(SP);
sz_t = size(t);

if sum(sz_SP - sz_t) ~= 0 
    error('eos80_legacy_neutral_surfaces: SP and t need be of the same size')
end 

p = gsw_resize(p,sz_SP);

sz_gn = size(gamma_n);

if sum(sz_SP - sz_gn) ~= 0 
    error('eos80_legacy_neutral_surfaces: SP and gamma_n need be of the same size')
end 

% Assign a bin for each bottle.
[profile_length, number_of_profiles] = size(SP);
number_of_surfaces = length(gamma_n_surfaces);
Ishallower = NaN(number_of_surfaces,number_of_profiles);

for I = 1:number_of_surfaces
    
    start_point = (I-1)*number_of_profiles + 1;
    end_point = (I)*number_of_profiles;
    
    gamma_diff = gamma_n - gamma_n_surfaces(I);

    gamma_diff_shallower = gamma_diff;
    gamma_diff_shallower(gamma_diff>0) = NaN;
    
    [dummy, Ibottle] = max(gamma_diff_shallower);
    Ibottle(isnan(dummy)) = NaN;
    
    Ishallower(I,:) = Ibottle(:).';
    Ishallow(start_point:end_point) = sub2ind(size(gamma_n),Ishallower(I,:),[1:number_of_profiles]);
    Ideep(start_point:end_point) = sub2ind(size(gamma_n),Ishallower(I,:)+1,[1:number_of_profiles]);
    
end
    
[Inn] = find(~isnan(Ishallow) & ~isnan(Ideep));
p_mid = 0.5.*(p(Ishallow(Inn)) + p(Ideep(Inn)));

r = (p_mid - p(Ishallow(Inn)))./(p(Ideep(Inn)) - p(Ishallow(Inn)));

SP_mid = SP(Ishallow(Inn)) + r.*(SP(Ideep(Inn)) - SP(Ishallow(Inn)));

pt0_shallow = eos80_legacy_pt(SP(Ishallow(Inn)),t(Ishallow(Inn)),p(Ishallow(Inn)),0);
pt0_deep = eos80_legacy_pt(SP(Ideep(Inn)),t(Ideep(Inn)),p(Ideep(Inn)),0);
pt0_mid = pt0_shallow + r.*(pt0_deep - pt0_shallow);
t_mid = eos80_legacy_pt(SP_mid,pt0_mid,0,p_mid);

[alpha_shallow, beta_shallow] = eos80_legacy_alpha_beta_wrt_pt(SP(Ishallow(Inn)),t(Ishallow(Inn)),p(Ishallow(Inn)));
[alpha_deep, beta_deep] = eos80_legacy_alpha_beta_wrt_pt(SP(Ideep(Inn)),t(Ideep(Inn)),p(Ideep(Inn)));

alpha_mid = -0.5.*(alpha_shallow + alpha_deep);
beta_mid = 0.5.*(beta_shallow + beta_deep);
[ma,na] = size(alpha_mid);

sigma_mid = eos80_legacy_sigma(SP_mid,t_mid,p_mid);
[mra,nra] = size(sigma_mid);
if mra == na & nra == ma
    sigma_mid = sigma_mid.';
end

rho_mid = 1000 + sigma_mid;

pt_shallow = eos80_legacy_pt(SP(Ishallow(Inn)),t(Ishallow(Inn)),p(Ishallow(Inn)),0);
pt_deep = eos80_legacy_pt(SP(Ideep(Inn)),t(Ideep(Inn)),p(Ideep(Inn)),0);

delta_SP = SP(Ideep(Inn)) - SP(Ishallow(Inn));
delta_pt = pt_deep - pt_shallow;

p_shallow = p(Ishallow(Inn));
p_deep = p(Ideep(Inn));

delta_p = p_deep - p_shallow;
delta_p2 = delta_p.*delta_p;
            
bden = rho_mid.*(beta_mid.*delta_SP - alpha_mid.*delta_pt);
bden(abs(bden) < 1e-6) = 1e-6;

b_mid = (gamma_n(Ideep(Inn)) - gamma_n(Ishallow(Inn)))./bden;

a = delta_SP.*(beta_deep - beta_shallow) - delta_pt.*((-1).*alpha_deep - (-1).*alpha_shallow);
a = (a.*b_mid.*rho_mid)./(2.*delta_p2);

b = delta_SP.*(p_deep.*beta_shallow - p_shallow.*beta_deep) - delta_pt.*(p_deep.*(-1).*alpha_shallow - p_shallow.*(-1).*alpha_deep);
b = (b.*b_mid.*rho_mid)./delta_p2;

c = NaN(1,number_of_surfaces.*number_of_profiles); 
c_part = delta_SP.*(beta_shallow.*(p_shallow - 2.*p_deep) + beta_deep.*p_shallow) - delta_pt.*((-1).*alpha_shallow.*(p_shallow - 2.*p_deep) + (-1).*alpha_deep.*p_shallow);
c(Inn) = gamma_n(Ishallow(Inn)) + (b_mid.*rho_mid.*p_shallow.*c_part)./(2.*delta_p2);
    
[mc,nc] = size(c);
if mc == na & nc == ma
    c = c.';
end

for I = 1:number_of_surfaces
    start_point = (I-1).*number_of_profiles + 1;
    end_point = I.*number_of_profiles;
    c(start_point:end_point) = c(start_point:end_point) - gamma_n_surfaces(I);
end

p_tol = 1e-3;

p_ns = NaN(1,number_of_surfaces.*number_of_profiles);
[mp_ns,np_ns] = size(p_ns);
if mp_ns == na & np_ns == ma
    p_ns = p_ns.';
end
SP_ns = p_ns;
t_ns = p_ns;

if any(a(:) ~=0 & bden(:) ~= 1e-6)
    [Ipart] = find(a ~=0 & bden ~= 1e-6);
    
    q_part1 = b(Ipart).*b(Ipart) - 4.*a(Ipart).*c(Inn(Ipart));
    q_part1(q_part1 < 0) = NaN;
    
    q = -(b(Ipart) + sign(b(Ipart)).*sqrt(q_part1)).*0.5;
    
    pns1 = q./a(Ipart);
    pns2 = c(Inn(Ipart))./q;
    
    if any(pns1 >= (p_shallow(Ipart) - p_tol) & pns1 <= (p_deep(Ipart) + p_tol))
        [Ipart2] = find(pns1 >= (p_shallow(Ipart) - p_tol) & pns1 <= (p_deep(Ipart) + p_tol));       
        p_ns(Inn(Ipart(Ipart2))) = min(p_deep(Ipart(Ipart2)),max(pns1(Ipart2),p_shallow(Ipart(Ipart2))));
    end
    
    if any(pns2 >= (p_shallow(Ipart) - p_tol) & pns2 <= (p_deep(Ipart) + p_tol)) 
        [Ipart2] = find(pns2 >= (p_shallow(Ipart)-p_tol) & pns2 <= (p_deep(Ipart) + p_tol));
        p_ns(Inn(Ipart(Ipart2))) = min(p_deep(Ipart(Ipart2)),max(pns2(Ipart2),p_shallow(Ipart(Ipart2))));
    end    
end

if any(a(:) ==0 & bden(:) == 1e-6)
    [Ipart] = find(a ==0 & bden == 1e-6);
    rg = NaN(size(p_ns));
    
    for I = 1:number_of_surfaces
        start_point = (I-1).*number_of_profiles + 1;
        end_point = I.*number_of_profiles;
        
        rg(start_point:end_point) = (gamma_n_surfaces(I) - gamma_n(Ishallow(start_point:end_point)))./(gamma_n(Ideep(start_point:end_point)) - gamma_n(Ishallow(start_point:end_point)));      
    end
    
    p_ns(Ipart) = p_shallow(Ipart) + rg(Inn(Ipart)).*(p_deep(Ipart) - p_shallow(Ipart));
end

r = (p_ns(Inn) - p_shallow)./(p_deep - p_shallow);

SP_ns(Inn) = SP(Ishallow(Inn)) + r.*(SP(Ideep(Inn)) - SP(Ishallow(Inn)));

pt0_shallow = eos80_legacy_pt(SP(Ishallow(Inn)),t(Ishallow(Inn)),p(Ishallow(Inn)),0);
pt0_deep = eos80_legacy_pt(SP(Ideep(Inn)),t(Ideep(Inn)),p(Ideep(Inn)),0);
pt0_ns = pt0_shallow + r.*(pt0_deep - pt0_shallow);
t_ns(Inn) = eos80_legacy_pt(SP_ns(Inn),pt0_ns,0,p_ns(Inn));

dims = length(size(SP));
if dims == 2
    SP_neutral_surface = NaN(number_of_surfaces,number_of_profiles);
else
    [dummy,number_of_profiles1,number_of_profiles2] = size(SP);
    SP_neutral_surface = NaN(number_of_surfaces,number_of_profiles1,number_of_profiles2);
end
t_neutral_surface = SP_neutral_surface;
p_neutral_surface = SP_neutral_surface;

for I = 1:number_of_surfaces
    start_point = (I-1).*number_of_profiles + 1;
    end_point = I.*number_of_profiles;
    SP_neutral_surface(I,:) = SP_ns(start_point:end_point);
    t_neutral_surface(I,:) = t_ns(start_point:end_point);
    p_neutral_surface(I,:) = p_ns(start_point:end_point);
end

end
