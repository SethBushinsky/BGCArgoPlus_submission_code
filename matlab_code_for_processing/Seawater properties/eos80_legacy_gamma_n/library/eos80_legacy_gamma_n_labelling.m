function [gamma_n, gamma_error_lower, gamma_error_upper] = eos80_legacy_gamma_n_labelling(SP,t,p,long,lat)

% eos80_legacy_gamma_n_labelling                            neutral density
%==========================================================================
%
% USAGE:
%  [gamma_n, gamma_error_lower, gamma_error_upper] = ...
%                    eos80_legacy_gamma_n_labelling(SP,t,p,long,lat)
%
% DESCRIPTION:
%  Calculates neutral density value, gamma_n, in the open ocean by 
%  spatially interpolating the global reference data set, gamma_n_ref, to
%  the location of the seawater sample using the neutral tangent plane.
%  
%  The correct way of calculating gamma_n is by calling 
%  eos80_legacy_gamma_n.  
%
%  This function uses version 1.0 of the gamma_n look up table (1997). 
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
%  SP, t, p, long & lat need to be vectors and have the same dimensions.
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
% VERSION NUMBER: 3.05.9 (1st March, 2017)
%
% REFERENCES:
%  Jackett, D.R., and T.J. McDougall, 1997: A Neutral Density Variable for
%   the World's Oceans. J. Phys. Oceanogr., 27, pp. 237–263.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Start of the calculation (extracting from a look up table)
%--------------------------------------------------------------------------
persistent gamma_n_ref lat_ref long_ref p_ref SP_ref t_ref a_ref

if isempty(gamma_n_ref)
    gamma_n_data = 'eos80_legacy_gamma_n_data.mat';
    
    gamma_n_data_file = which(gamma_n_data);

    load (gamma_n_data_file,'gamma_n_ref','lat_ref','long_ref','p_ref',...
        'a_ref','SP_ref','t_ref','n_ref');
end

% precalculate constants 
nx = length(long_ref); 
ny = length(lat_ref); 
nz = length(p_ref); 

% Calculate the 4 grid points surrounding the data
Ix0 = floor(1 + (nx-1)*(long - long_ref(1))./(long_ref(nx) - long_ref(1)));
Ix0 = Ix0(:); 
Ix0(Ix0 == nx) = nx - 1;
              
Iy0 = floor(1 + (ny-1)*(lat - lat_ref(1))./(lat_ref(ny) - lat_ref(1)));
Iy0 = Iy0(:); 
Iy0(Iy0 == ny) = ny - 1;
     
Iy0_Ix0ny = Iy0 + Ix0.*ny;        
I1 = Iy0_Ix0ny - ny;              
I2 = Iy0_Ix0ny;
I3 = Iy0_Ix0ny + 1;
I4 = Iy0_Ix0ny + (1 - ny);

% calculate the ntp for the bottle to each of the 4 surrounding grid points
[mp,np] = size(I1);
p_ref1_dummy = repmat(p_ref(:),1,mp);
p_ref1_dummy(isnan(SP_ref(:,I1))) = NaN;
[mp,np] = size(I2);
p_ref2_dummy = repmat(p_ref(:),1,mp);
p_ref2_dummy(isnan(SP_ref(:,I2))) = NaN;
[mp,np] = size(I3);
p_ref3_dummy = repmat(p_ref(:),1,mp);
p_ref3_dummy(isnan(SP_ref(:,I3))) = NaN;
[mp,np] = size(I4);
p_ref4_dummy = repmat(p_ref(:),1,mp);
p_ref4_dummy(isnan(SP_ref(:,I4))) = NaN;

% do the quick vectorised ntp calculation
[SP_ntp1, t_ntp1, p_ntp1] = eos80_legacy_ntp_bottle_to_cast(SP(:).',t(:).',p(:).',SP_ref(:,I1),t_ref(:,I1),p_ref1_dummy);
[SP_ntp2, t_ntp2, p_ntp2] = eos80_legacy_ntp_bottle_to_cast(SP(:).',t(:).',p(:).',SP_ref(:,I2),t_ref(:,I2),p_ref2_dummy);
[SP_ntp3, t_ntp3, p_ntp3] = eos80_legacy_ntp_bottle_to_cast(SP(:).',t(:).',p(:).',SP_ref(:,I3),t_ref(:,I3),p_ref3_dummy);
[SP_ntp4, t_ntp4, p_ntp4] = eos80_legacy_ntp_bottle_to_cast(SP(:).',t(:).',p(:).',SP_ref(:,I4),t_ref(:,I4),p_ref4_dummy);

% Assign a pressure bin for each bottle.
n1 = length(p_ntp1);
Iz1 = NaN(n1,1);
Iz2 = Iz1;
Iz3 = Iz1;
Iz4 = Iz1;
for I = 2:nz   
    Iz1(p_ntp1 >= p_ref(I-1) & p_ntp1 < p_ref(I)) = I - 1;    
    Iz2(p_ntp2 >= p_ref(I-1) & p_ntp2 < p_ref(I)) = I - 1;    
    Iz3(p_ntp3 >= p_ref(I-1) & p_ntp3 < p_ref(I)) = I - 1;    
    Iz4(p_ntp4 >= p_ref(I-1) & p_ntp4 < p_ref(I)) = I - 1;     
end
Iz1(p_ntp1 == p_ref(nz)) = nz; 
Iz2(p_ntp2 == p_ref(nz)) = nz; 
Iz3(p_ntp3 == p_ref(nz)) = nz; 
Iz4(p_ntp4 == p_ref(nz)) = nz; 

% calculate the indicies for each of the 4 surrounding non-NaN npt's. 
I1z_nn = find(~isnan(Iz1));
Ixyz1 = sub2ind(size(SP_ref),Iz1(I1z_nn),Iy0(I1z_nn),Ix0(I1z_nn));

I2z_nn = find(~isnan(Iz2));
Ixyz2 = sub2ind(size(SP_ref),Iz2(I2z_nn),Iy0(I2z_nn),Ix0(I2z_nn)+1);

I3z_nn = find(~isnan(Iz3));
Ixyz3 = sub2ind(size(SP_ref),Iz3(I3z_nn),Iy0(I3z_nn)+1,Ix0(I3z_nn)+1);

I4z_nn = find(~isnan(Iz4));
Ixyz4 = sub2ind(size(SP_ref),Iz4(I4z_nn),Iy0(I4z_nn)+1,Ix0(I4z_nn));

% calculate a gamma for the bottle to each of the 4 surrounding grid points
gamma_ntp1 = NaN(1,length(SP));
rec_dp = 1./(p_ref(Iz1(I1z_nn)+1) - p_ref(Iz1(I1z_nn)));
p1 = ((p_ntp1(I1z_nn) - p_ref(Iz1(I1z_nn))).*rec_dp).';
p2 = ((p_ntp1(I1z_nn) - p_ref(Iz1(I1z_nn)+1)).*rec_dp).';
gamma_ntp1(I1z_nn) = gamma_n_ref(Ixyz1) + p1.*(gamma_n_ref(Ixyz1+1) - gamma_n_ref(Ixyz1) + p2.*a_ref(Ixyz1));

gamma_ntp2 = NaN(1,length(SP));
rec_dp = 1./(p_ref(Iz2(I2z_nn)+1) - p_ref(Iz2(I2z_nn)));
p1 = ((p_ntp2(I2z_nn) - p_ref(Iz2(I2z_nn))).*rec_dp).';
p2 = ((p_ntp2(I2z_nn) - p_ref(Iz2(I2z_nn)+1)).*rec_dp).';
gamma_ntp2(I2z_nn) = gamma_n_ref(Ixyz2) + p1.*(gamma_n_ref(Ixyz2+1)  - gamma_n_ref(Ixyz2) + p2.*a_ref(Ixyz2));

gamma_ntp3 = NaN(1,length(SP));
rec_dp = 1./(p_ref(Iz3(I3z_nn)+1) - p_ref(Iz3(I3z_nn)));
p1 = ((p_ntp3(I3z_nn) - p_ref(Iz3(I3z_nn))).*rec_dp).';
p2 = ((p_ntp3(I3z_nn) - p_ref(Iz3(I3z_nn)+1)).*rec_dp).';
gamma_ntp3(I3z_nn) = gamma_n_ref(Ixyz3) + p1.*(gamma_n_ref(Ixyz3+1) - gamma_n_ref(Ixyz3) + p2.*a_ref(Ixyz3));

gamma_ntp4 = NaN(1,length(SP));
rec_dp = 1./(p_ref(Iz4(I4z_nn)+1) - p_ref(Iz4(I4z_nn)));
p1 = ((p_ntp4(I4z_nn) - p_ref(Iz4(I4z_nn))).*rec_dp).';
p2 = ((p_ntp4(I4z_nn) - p_ref(Iz4(I4z_nn)+1)).*rec_dp).';
gamma_ntp4(I4z_nn) = gamma_n_ref(Ixyz4) + p1.*(gamma_n_ref(Ixyz4+1) - gamma_n_ref(Ixyz4) + p2.*a_ref(Ixyz4));

if nargout > 1
    % calculate the errors, when there is a gamma value
    g1_error = nan(4,length(gamma_ntp1));
    g2_lower_error = g1_error;
    g2_upper_error = g1_error;
    if any(~isnan(gamma_ntp1))
        [Inn] = find(~isnan(gamma_ntp1));
        [g1_error(1,Inn),g2_lower_error(1,Inn),g2_upper_error(1,Inn)] = ...
            eos80_legacy_gamma_n_labelling_error(SP(Inn).',t(Inn).',p(Inn).',lat(Inn).',...
            SP_ref(:,I1(Inn)),t_ref(:,I1(Inn)),p_ref1_dummy(:,Inn),gamma_n_ref(:,I1(Inn)),a_ref(:,I1(Inn)),...
            SP_ntp1(Inn),t_ntp1(Inn),p_ntp1(Inn),gamma_ntp1(Inn));
    end
    if any(~isnan(gamma_ntp2))
        [Inn] = find(~isnan(gamma_ntp2));
        [g1_error(2,Inn),g2_lower_error(2,Inn),g2_upper_error(2,Inn)] = ...
            eos80_legacy_gamma_n_labelling_error(SP(Inn).',t(Inn).',p(Inn).',lat(Inn).',...
            SP_ref(:,I2(Inn)),t_ref(:,I2(Inn)),p_ref2_dummy(:,Inn),gamma_n_ref(:,I2(Inn)),a_ref(:,I2(Inn)),...
            SP_ntp2(Inn),t_ntp2(Inn),p_ntp2(Inn),gamma_ntp2(Inn));
    end
    if any(~isnan(gamma_ntp3))
        [Inn] = find(~isnan(gamma_ntp3));
        [g1_error(3,Inn),g2_lower_error(3,Inn),g2_upper_error(3,Inn)] = ...
            eos80_legacy_gamma_n_labelling_error(SP(Inn).',t(Inn).',p(Inn).',lat(Inn).',...
            SP_ref(:,I3(Inn)),t_ref(:,I1(Inn)),p_ref3_dummy(:,Inn),gamma_n_ref(:,I3(Inn)),a_ref(:,I3(Inn)),...
            SP_ntp3(Inn),t_ntp3(Inn),p_ntp3(Inn),gamma_ntp3(Inn));
    end
    if any(~isnan(gamma_ntp4))
        [Inn] = find(~isnan(gamma_ntp4));
        [g1_error(4,Inn),g2_lower_error(4,Inn),g2_upper_error(4,Inn)] = ...
            eos80_legacy_gamma_n_labelling_error(SP(Inn).',t(Inn).',p(Inn).',lat(Inn).',...
            SP_ref(:,I4(Inn)),t_ref(:,I4(Inn)),p_ref4_dummy(:,Inn),gamma_n_ref(:,I4(Inn)),a_ref(:,I4(Inn)),...
            SP_ntp4(Inn),t_ntp4(Inn),p_ntp4(Inn),gamma_ntp4(Inn));
    end
end

extrap_weights = ones(4,length(gamma_ntp1));

%calculate a gamma and the appropiate weights, when extrapolation is needed
if any(isnan(gamma_ntp1))
    [Inan] = find(isnan(gamma_ntp1));
    [gamma_ntp1(Inan),extrap_weights(1,Inan)] = ...
        eos80_legacy_gamma_n_labelling_extrapolation(SP(Inan),t(Inan),p(Inan),...
        SP_ref(:,I1(Inan)),t_ref(:,I1(Inan)),p_ref(:),gamma_n_ref(:,I1(Inan)));
end
if any(isnan(gamma_ntp2))
    [Inan] = find(isnan(gamma_ntp2));
    [gamma_ntp2(Inan), extrap_weights(2,Inan)] = ...
        eos80_legacy_gamma_n_labelling_extrapolation(SP(Inan),t(Inan),p(Inan),...
        SP_ref(:,I2(Inan)),t_ref(:,I2(Inan)),p_ref(:),gamma_n_ref(:,I2(Inan)));
end
if any(isnan(gamma_ntp3))
    [Inan] = find(isnan(gamma_ntp3));
    [gamma_ntp3(Inan), extrap_weights(3,Inan)] = ...
        eos80_legacy_gamma_n_labelling_extrapolation(SP(Inan),t(Inan),p(Inan),...
        SP_ref(:,I3(Inan)),t_ref(:,I3(Inan)),p_ref(:),gamma_n_ref(:,I3(Inan)));
end
if any(isnan(gamma_ntp4))
    [Inan] = find(isnan(gamma_ntp4));
    [gamma_ntp4(Inan), extrap_weights(4,Inan)] = ...
        eos80_legacy_gamma_n_labelling_extrapolation(SP(Inan),t(Inan),p(Inan),...
        SP_ref(:,I4(Inan)),t_ref(:,I4(Inan)),p_ref(:),gamma_n_ref(:,I4(Inan)));
end

gamma_ntps(1,:) = gamma_ntp1;
gamma_ntps(2,:) = gamma_ntp2;
gamma_ntps(3,:) = gamma_ntp3;
gamma_ntps(4,:) = gamma_ntp4;

% calculate spatial weights, based on their distance from each grid point
rx = (long(:).' - long_ref(Ix0))./(long_ref(Ix0+1) - long_ref(Ix0));
ry = (lat(:).' - lat_ref(Iy0))./(lat_ref(Iy0+1) - lat_ref(Iy0));

% separate the oceans either side of Central America
if any(abs(long-277.6085)<=17.6085 & abs(lat-9.775) <= 9.775)
    [I_ca] = find(abs(long-277.6085)<=17.6085 & abs(lat-9.775) <= 9.775);
    gamma_ntps(:,I_ca) = gsw_refdata_add_ca_barrier(gamma_ntps(:,I_ca),  long(I_ca), lat(I_ca),...
             long_ref(Ix0(I_ca)), lat_ref(Iy0(I_ca)), 4, 4, extrap_weights(:,I_ca));
end

% replace any NaN's with the mean of the non-NaN grid points
if any(isnan(sum(gamma_ntps))) | any(sum(extrap_weights) < 4)
    [I_nan] = find(isnan(sum(gamma_ntps)) | sum(extrap_weights) < 1);
    gamma_ntps(:,I_nan) = gsw_refdata_replace_nan(gamma_ntps(:,I_nan), extrap_weights(:,I_nan));
end

% calculate the averaged gamma_n
gamma_n = (1-ry).*(gamma_ntps(1,:) + rx.*(gamma_ntps(2,:) - gamma_ntps(1,:))) ...
    + ry.*(gamma_ntps(4,:) + rx.*(gamma_ntps(3,:) - gamma_ntps(4,:)));

if nargout > 1
    % calculate the averaged gamma_n_error
    gamma_error1 = (1-ry).*(g1_error(1,:) + rx.*(g1_error(2,:) - g1_error(1,:))) ...
        + ry.*(g1_error(4,:) + rx.*(g1_error(3,:) - g1_error(4,:)));
    
    % calculate the averaged difference between the averaged gamma_n and the individual gamma_n's
    gamma_error3 = (1-ry).*(abs(gamma_ntps(1,:) - gamma_n) ...
        + rx.*(abs(gamma_ntps(2,:) - gamma_n) - abs(gamma_ntps(1,:) - gamma_n) )) ...
        + ry.*(abs(gamma_ntps(4,:) - gamma_n) ...
        + rx.*(abs(gamma_ntps(3,:) - gamma_n) - abs(gamma_ntps(4,:) - gamma_n)));
    
    gamma_errors(1,:) = gamma_error1;
    gamma_errors(2:5,:) = g2_lower_error;
    gamma_errors(6,:) = gamma_error3;
    
    gamma_error_lower = max(gamma_errors);
    gamma_error_lower = gamma_error_lower(:);
    
    gamma_errors(2:5,:) = g2_upper_error;
    gamma_error_upper = max(gamma_errors);
    gamma_error_upper = gamma_error_upper(:);
end

gamma_n = gamma_n(:);

end
