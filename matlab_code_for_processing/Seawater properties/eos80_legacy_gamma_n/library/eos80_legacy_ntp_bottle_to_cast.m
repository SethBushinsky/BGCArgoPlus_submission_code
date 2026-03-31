function [SP_ntp, t_ntp, p_ntp, p_bar] = eos80_legacy_ntp_bottle_to_cast(SP_bottle,t_bottle,p_bottle,SP,t,p)

% eos80_legacy_ntp_bottle_to_cast               Practical Salinity, in-situ
%                     temperature and pressure on a neighbouring cast using
%                                the neutral tangent plane (eos80 equation)
%==========================================================================
%
% USAGE:
%  [SP_ntp, t_ntp, p_ntp, p_bar] = ...
%       gsw_ntp_bottle_to_cast(SP_bottle,t_bottle,p_bottle,SA,CT,p)
%
% DESCRIPTION:
%  This function finds the Practical Salinity, in-situ temperature and
%  pressure on the casts, which have the same potential density referred to
%  the mid-point pressure, p_bar, as the bottle.
%
%  This programme assumes that the neighbouring cast is stable.  It
%  searches in one direction only, based on the initial value of the
%  difference in the locally referenced potential density.  It searches for
%  the corresponding locally referenced potential density on the 
%  neighbouring cast starting from the bottle pressure.
%
%  Note that this function is intended to work on stable profiles.  If
%  there are instabillities then there are mupliple solutions. It will
%  not find all the crossings. The results may including stable, unstable
%  and neutral crossings.
%
% INPUT:
%  SP_bottle =  Practical Salinity of the bottle               [ unitless ]
%  t_bottle  =  in-situ temperature of the bottle(ITS-90)         [ deg C ]
%  p_bottle  =  sea pressure of the bottle                         [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  SP        =  Practical Salinity of the neighbouring cast    [ unitless ]
%  t         =  in-situ temperature of the neighbouring cast      [ deg C ]
%  p         =  sea pressure of the neighbouring cast              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SP_bottle, t_bottle and p_bottle need to be single values, having
%  dimensions 1x1.
%  SP, t and p need to have the same dimensions Mx1 or MxN, where M is the
%  number of bottles in the cast and N is the number of profiles.
%  Note that they need to be arranged in a column, if they are inputed as a
%  row the programme will assume they are seperate casts.
%
% OUTPUT:
%  SP_ntp = Practical Salinity of the neutral tangent plane    [ unitless ]
%  t_ntp  = in-situ Temperature of the neutral tangent plane      [ deg C ]
%  p_ntp  = pressure of the neutral tangent plane                  [ dbar ]
%  p_bar  = the average pressure of the neutral tangent plane      [ dbar ]
%
% AUTHOR: David Jackett
%  Modified by Guillaume Serazin, Stefan Riha and Paul Barker
%
% VERSION NUMBER: 3.06 (18th September, 2017)
%
% REFERENCES:
%  Jackett, D.R., and T.J. McDougall, 1997: A Neutral Density Variable for
%   the World's Oceans. J. Phys. Oceanogr., 27, pp. 237–263.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

delta = 1e-9; % The accepted tolerance when calculating the locally
% referenced potential density difference between the bottle and the
% neighbouring cast.  Note that this is approximately equal to 0.1 dbar.

[profile_length, Number_of_profiles] = size(SP);
size_dummy = size(SP);

Idata = 1:Number_of_profiles;

size_dummy(1) = 1;
SP_ntp = NaN(size_dummy);
t_ntp = SP_ntp;
p_ntp = SP_ntp;
p_bar = SP_ntp;

% discard NaN's (land), Idata may contain the indices of remaining points
[SP_bottle, t_bottle, p_bottle, SP, t, p, Idata] = ...
    gsw_gamma_discard_NaNs(SP_bottle, t_bottle, p_bottle, SP, t, p, Idata);

[profile_length, Number_of_profiles] = size(SP);

Iwet = [1:Number_of_profiles];

p0_stacked = repmat(p_bottle(:)',[profile_length 1]);

I_bottle_above = sum(p0_stacked >= p,1);                         % Index of adjacent bottle looking up, values 1 to nz.
bottles_above = p_bottle(:)' < p(1,:);                           % bottle above shallowest cast data point
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
    error('eos80_legacy_ntp_bottle_to_cast: error # 1. Email help@teos-10.org')          % land should have been discarded above
end

SP_dummy = SP(end:-1:1,:,:);
Idata_nn = ~isnan(SP_dummy);
Idata_intergral = cumsum(Idata_nn,1);
su = ones(size(Idata_intergral));
su(Idata_intergral~=0) = 0;
Ideepest = sum(su,1) + 1;
Ideepest = profile_length - Ideepest + 1;        % vertical index of deepest data point
if any(Ideepest < 1 | Ideepest > profile_length)
    error('eos80_legacy_ntp_bottle_to_cast: error # 2. Email help@teos-10.org')           % land should have been discarded above
end

go_shallow_old = false(1,Number_of_profiles);
go_deep_old = false(1,Number_of_profiles);

bisect = false(1,Number_of_profiles);
done = false(1,Number_of_profiles);
Number_of_iterations = 0;

while any(~done)
    
    Number_of_iterations = Number_of_iterations + 1;
    
    p_mid = 0.5*(p_bottle + p_tmp);
        
    pot_rho_bottle = eos80_legacy_pot_rho_pref(SP_bottle,t_bottle,p_bottle,p_mid);
    pot_rho_cast = eos80_legacy_pot_rho_pref(SP_tmp,t_tmp,p_tmp,p_mid);
    pot_rho_diff =  pot_rho_cast - pot_rho_bottle;
    
    go_shallow = (pot_rho_diff >=  delta);   % The intercept is shallower
    go_deep = (pot_rho_diff <= -delta);      % The intercept is deeper
        
    searched(~isnan(pot_rho_diff)) = true;
    
    % If pot_rho_diff, the locally referenced potential density difference between
    % the bottle and the neighbouring cast, is NaN at the current pressure
    % but it is well defined elsewhere in the watercolumn.
    nskip_up = I_bottle_above(Iwet) - Ideepest(Iwet);
    nskip_down = Ishallowest(Iwet) - I_bottle_above(Iwet);
    
    go_shallow_nan = isnan(pot_rho_diff)  & ~searched & (nskip_up>0);
    go_deep_nan = isnan(pot_rho_diff) & ~searched & (nskip_down>0);
    
    go_shallow = go_shallow | go_shallow_nan;
    go_deep = go_deep | go_deep_nan;
    
    final = abs(pot_rho_diff) < delta;  % The calculated salinty, temperature and
    % pressure are within the tolerence level.
    
    SP_ntp(Idata(Iwet(final))) = SP_tmp(final);
    t_ntp(Idata(Iwet(final))) = t_tmp(final);
    p_ntp(Idata(Iwet(final))) = p_tmp(final);
    p_bar(Idata(Iwet(final))) = p_mid(final);
    
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
    
    done = (isnan(pot_rho_diff) & searched) | final | out2;
    
    I3d_wet = I3d_wet(~out);
    search_initial = search_initial(~done);
    bisect = bisect(~done);
    crossed = crossed(~done);
    searched = searched(~done);
    
    if ~all((search_initial & ~bisect) | (~search_initial & bisect)) % either bisect or keep searching
        error('eos80_legacy_ntp_bottle_to_cast: # 3. Email help@teos-10.org)')
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
    % potential density between the bottle and the cast
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
