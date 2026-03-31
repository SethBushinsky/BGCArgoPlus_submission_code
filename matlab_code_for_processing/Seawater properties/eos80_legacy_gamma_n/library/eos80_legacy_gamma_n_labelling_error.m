function [g1_error,g2_lower_error,g2_upper_error] = eos80_legacy_gamma_n_labelling_error(SP,t,p,lat,SP_ref,t_ref,p_ref,gamma_n_ref,a_ref,SP_ntp,t_ntp,p_ntp,gamma_ntp)

% eos80_legacy_gamma_n_labelling_error     neutral density labelling errors
%==========================================================================
%
% USAGE:
%  [g1_error, g2_lower_error, g2_upper_error] = ...
%    eos80_legacy_gamma_n_labelling_error(SP,t,p,lat,...
%               SP_ref,t_ref,p_ref,gamma_n_ref,a_ref,...
%               SP_ntp,t_ntp,p_ntp,gamma_ntp)
%
% DESCRIPTION:
%  Calculates an estimate of the neutral density errors.
%
% INPUT:
%  SP        = Practical Salinity of the bottle                [ unitless ]
%  t         = in-situ temperature (ITS-90) of the bottle         [ deg C ]
%  p         = sea pressure of the bottle                          [ dbar ] 
%          ( i.e. absolute pressure - 10.1325 dbar )
%  lat       = latitude in decimal degrees north            [ -90 ... +90 ]
%  SP_ref    = Practical Salinity of the reference cast        [ unitless ]
%  t_ref     = in-situ temperature of the reference cast          [ deg C ]
%  p_ref     = sea pressure of the reference cast                  [ dbar ]
%  gamma_n_ref = neutral density of the reference cast           [ kg/m^3 ]
%  a_ref     = quadratic coefficient of the reference cast
%  SP_ntp    = Practical Salinity on the neutral tangent plane [ unitless ]
%  t_ntp     = in-situ temperature on the neutral tangent plane   [ deg C ]
%  p_ntp     = sea pressure on the neutral tangent plane           [ dbar ]
%  gamma_ntp = neutral density on the neutral tangent plane      [ kg/m^3 ]
%
%  SP, t, p, long need to be vectors and have the same dimensions.
%
% OUTPUT:
%  g1_error        =  p-theta error                                [ dbar ]
%  g2_lower_error  =  lower submesoscale coherent vortices gamma error
%                                                                  [ dbar ]
%  g2_upper_error  =  upper submesoscale coherent vortices gamma error
%                                                                  [ dbar ]
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
%  McDougall, T.J., 1987: The vertical motion of submesoscale coherent
%   vortices across neutral surfaces. Journal of physical oceanography, 
%   17, pp.2334-2342.
%
% Unesco, 1983: Algorithms for computation of fundamental properties of 
%  seawater. Unesco technical papers in marine science 44, 53pp.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

p0 = 0;
Tb = 2.7e-8;
gamma_limit = 26.845;
test_limit = 0.1;

% p-theta error
pt = eos80_legacy_pt(SP,t,p,p0);
pt_ntp = eos80_legacy_pt(SP_ntp,t_ntp,p_ntp,p0);

rho_ntp = 1000 + eos80_legacy_sigma(SP_ntp,t_ntp,p_ntp);

Ins = NaN(1,length(gamma_ntp));
for I = 1:length(gamma_ntp)
    [dummy, Ins(I)] = max(gamma_n_ref(gamma_n_ref(:,I) < gamma_ntp(I)));
end

pmid = 0.5.*(p_ref(Ins) + p_ref(Ins+1));
sig_l = eos80_legacy_pot_rho_pref(SP_ref(Ins),t_ref(Ins),p_ref(Ins),pmid);
sig_h = eos80_legacy_pot_rho_pref(SP_ref(Ins+1),t_ref(Ins+1),p_ref(Ins+1),pmid);

b = (gamma_n_ref(Ins+1) - gamma_n_ref(Ins))./(sig_h - sig_l);

dp = p_ntp - p;

delta_pt = pt_ntp - pt;

g1_error = rho_ntp.*b.*Tb.*abs(dp.*delta_pt)./6;

g2_lower_error = zeros(size(SP));
g2_upper_error = g2_lower_error;

if any(lat <= -60 | gamma_n_ref(1) >= gamma_limit)
    
    [I] = find(lat <= -60 | gamma_n_ref(1,:) >= gamma_limit);

    drldp = (sig_h(I) - sig_l(I))./(rho_ntp(I).*(p_ref(Ins(I)+1) - p_ref(Ins(I))));
    
    test = Tb.*delta_pt(I)./drldp;

    if any(abs(test) >= test_limit)
        
        [I2] = find(abs(test) >= test_limit);
        pt_change = dp(I(I2)).*delta_pt(I(I2));
        
        if any(pt_change >= 0)
            [I3] = find(pt_change >= 0);   
            g2_upper_error(I(I2(I3))) = (3.*g1_error(I(I2(I3))))./(1 - test(I2(I3)));
        end
        
        if any(pt_change < 0)
            [I3] = find(pt_change < 0);
            g2_lower_error(I(I2(I3))) = (3.*g1_error(I(I2(I3))))./(1 - test(I2(I3)));
        end
    end 
    
    if any(abs(test) < test_limit)
        
        [I2] = find(abs(test) < test_limit);

        p_error_estimate = eos80_legacy_gamma_n_labelling_error_estimate(SP(I(I2)),t(I(I2)),p(I(I2)),SP_ref(:,I(I2)),t_ref(:,I(I2)),p_ref(:,I(I2)));

        % Assign a pressure bin for each bottle.
        [Inn] = find(~isnan(p_error_estimate));        
        Iscv = NaN(1,length(Inn));
        for I3 = 1:length(Inn)
            [dummy, Iscv(I3)] = max(p_ref(p_ref(:,I(I2(Inn(I3)))) < p_error_estimate(Inn(I3))));
        end

        [Iscv_shallow] = sub2ind(size(p_ref),Iscv,I(I2(Inn)));
        Iscv_deep = Iscv_shallow + 1;

        rec_dp = 1./(p_ref(Iscv_deep) - p_ref(Iscv_shallow));
        p1 = (((p_error_estimate(Inn)) - p_ref(Iscv_shallow)).*rec_dp);
        p2 = (((p_error_estimate(Inn)) - p_ref(Iscv_deep)).*rec_dp);
        gamma_scv = gamma_n_ref(Iscv_shallow) + p1.*(gamma_n_ref(Iscv_deep) - gamma_n_ref(Iscv_shallow) + p2.*a_ref(Iscv_shallow));

        if any(p_error_estimate(Inn) <= p_ntp(I(I2(Inn))))
            [I3] = find(p_error_estimate(Inn) <= p_ntp(I(I2(Inn))));
            g2_lower_error(I(I2(Inn(I3)))) = gamma_ntp(I(I2(Inn(I3)))) - gamma_scv(I3);
        end
        
        if any(p_error_estimate(Inn) > p_ntp(I(I2(Inn))))
            [I3] = find(p_error_estimate(Inn) > p_ntp(I(I2(Inn))));
            g2_upper_error(I(I2(Inn(I3)))) = gamma_scv(I3) - gamma_ntp(I(I2(Inn(I3))));
        end
        
    end
    
end
    
end

%##########################################################################

function p_scv = eos80_legacy_gamma_n_labelling_error_estimate(SP_bottle,t_bottle,p_bottle,SP,t,p)

delta = 1e-12; % The accepted tolerance when calculating the locally 
% referenced specific volume difference between the bottle and the 
% neighbouring cast.  Note that this is approximately equal to 0.1 dbar.

[profile_length, Number_of_profiles] = size(SP);
size_dummy = size(SP);

Idata = 1:Number_of_profiles;

size_dummy(1) = 1;
p_scv = NaN(size_dummy);

% discard NaN's (land), Idata may contain the indices of remaining points
[SP_bottle, t_bottle, p_bottle, SP, t, p, Idata] = gsw_gamma_discard_NaNs(SP_bottle, t_bottle, p_bottle, SP, t, p, Idata);

[profile_length, Number_of_profiles] = size(SP);

Iwet = [1:Number_of_profiles];

p0_stacked = repmat(p_bottle(:)',[profile_length 1]);

I_bottle_above = sum(p0_stacked >= p,1);                         % index of adjacent bottle looking up, values 1 to nz.
bottles_above = p_bottle(:)' < p(1,:);                           % bottle above shallowest cast data doint
I_bottle_above(bottles_above) = 1;                               % start at the uppermost cast pair
I3d = I_bottle_above + profile_length*(0:Number_of_profiles-1);  % changing to 3-d index, based on casts supplied.
Iinc = NaN(1,Number_of_profiles);                                % increment of the 3-d index

SP_tmp = SP(I3d);
t_tmp = t(I3d);
p_tmp = p(I3d);

searched = ~isnan(SP_tmp);

Idata_nn = ~isnan(SP); 
Idata_intergral = cumsum(Idata_nn,1);
su = ones(size(Idata_intergral));
su(Idata_intergral~=0) = 0;
Ishallowest = sum(su,1) + 1;                    % vertical index of shallowest bottle
if any(Ishallowest < 1 | Ishallowest > profile_length)
    error('eos80_legacy_gamma_n_labelling_error_estimate: error # 1. Email help@teos-10.org ')  % NaN's should have been discarded above
end

SP_dummy = SP(end:-1:1,:,:);
Idata_nn = ~isnan(SP_dummy);
Idata_intergral = cumsum(Idata_nn,1);
su = ones(size(Idata_intergral));
su(Idata_intergral~=0) = 0;
%su = zeros(size(Idata_intergral));
Ideepest = sum(su,1) + 1;
Ideepest = profile_length - Ideepest + 1;        % vertical index of deepest data point
if any(Ideepest < 1 | Ideepest > profile_length)
    error('eos80_legacy_gamma_n_labelling_error_estimate: error # 2. Email help@teos-10.org')    % NaN's should have been discarded above
end

go_shallow_old = false(1,Number_of_profiles);
go_deep_old = false(1,Number_of_profiles);

bisect = false(1,Number_of_profiles);
done = false(1,Number_of_profiles);
Number_of_iterations = 0;

while any(~done)
    
    Number_of_iterations = Number_of_iterations + 1;
    
    v_bottle = 1./eos80_legacy_sigma(SP_bottle,t_bottle,p_tmp);
    v_cast = 1./eos80_legacy_sigma(SP_tmp,t_tmp,p_tmp);  
    
    v_local_diff =  v_bottle - v_cast;

    go_shallow = (v_local_diff >=  delta);   % The intercept is shallower
    go_deep = (v_local_diff <= -delta);      % The intercept is deeper 
        
    searched(~isnan(v_local_diff)) = true;
    
    % If v_local_diff, the locally referenced specific volume difference between
    % the bottle and the neighbouring cast, is NaN at the current pressure 
    % but it is well defined elsewhere in the water column.
    nskip_up = I_bottle_above(Iwet) - Ideepest(Iwet);
    nskip_down = Ishallowest(Iwet) - I_bottle_above(Iwet);
    
    go_shallow_nan = isnan(v_local_diff)  & ~searched & (nskip_up>0);
    go_deep_nan = isnan(v_local_diff) & ~searched & (nskip_down>0);
    
    go_shallow = go_shallow | go_shallow_nan;
    go_deep = go_deep | go_deep_nan;
    
    final = abs(v_local_diff) < delta;  % The calculated salinty, temperature and 
                                        % pressure are within the tolerence level. 
    
    p_scv(Idata(Iwet(final))) = p_tmp(final);

    cfb = (go_shallow_old & go_deep & searched); % crossed from below
    cfa = (go_deep_old & go_shallow & searched); % crossed from above
    crossed = cfb | cfa;
    start_bis = (crossed & ~bisect);   % start bisection here
    bisect = (bisect | start_bis);     % bisect here
    
    search_initial = (go_deep | go_shallow) & ~bisect & ~final;
    
    Iinc(Iwet(go_deep & search_initial)) = 1;
    Iinc(Iwet(go_shallow & search_initial)) = -1;
    
    Iinc(Iwet(go_deep_nan & search_initial)) = nskip_down(Iwet(go_deep_nan & search_initial));
    Iinc(Iwet(go_shallow_nan & search_initial)) = -nskip_up(Iwet(go_shallow_nan & search_initial));
        
    I_bottle_above = I_bottle_above + Iinc;
    I3d = I3d + Iinc;
    Iwet_tmp = I_bottle_above(Iwet(search_initial));  
    I3d_wet = I3d(Iwet(search_initial)); 
    
    out = (Iwet_tmp < 1) | (Iwet_tmp > profile_length);   % above or below domain 
    out2 = false(1,length(final)); 
    out2(search_initial) = out; 
    
    done = (isnan(v_local_diff) & searched) | final | out2; 
    
    I3d_wet = I3d_wet(~out);
    search_initial = search_initial(~done);
    bisect = bisect(~done);
    crossed = crossed(~done);
    searched = searched(~done);
    
    if ~all((search_initial & ~bisect) | (~search_initial & bisect)) % either bisect or keep searching
        error('eos80_legacy_gamma_n_labelling_error_estimate: error # 3. Email help@teos-10.org)')
    end
    
    if Number_of_iterations > 1
        SP1 = SP1(~done);
        t1 = t1(~done);
        p1 = p1(~done);
        
        SP2 = SP2(~done);
        t2 = t2(~done);
        p2 = p2(~done);
        
        SP2(crossed) = SP1(crossed);
        t2(crossed) = t1(crossed);
        p2(crossed) = p1(crossed);
    else
        SP2 = SP_tmp(~done);
        t2 = t_tmp(~done);
        p2 = p_tmp(~done);
    end
    
    SP1 = SP_tmp(~done);
    t1 = t_tmp(~done);
    p1 = p_tmp(~done);
    
    % data for next evaluation of the difference in the locally referenced
    % specific volume between the bottle and the cast
    SP_bottle = SP_bottle(~done);
    t_bottle = t_bottle(~done);
    p_bottle = p_bottle(~done);
    
    SP_tmp = NaN(1,sum(~done));
    t_tmp = NaN(1,sum(~done));
    p_tmp = NaN(1,sum(~done));
    
    SP_tmp(search_initial) = SP(I3d_wet);
    t_tmp(search_initial) = t(I3d_wet);
    p_tmp(search_initial) = p(I3d_wet);
    
    if any(bisect)
        SP_tmp(bisect) = 0.5*(SP1(bisect) + SP2(bisect));
        t_tmp(bisect) = 0.5*(t1(bisect) + t2(bisect));
        p_tmp(bisect) = 0.5*(p1(bisect) + p2(bisect));
    end
    
    go_shallow_old = go_shallow(~done);
    go_deep_old = go_deep(~done);
    
    Iwet = Iwet(~done);
end

end

