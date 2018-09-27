function inches_rainfall = rain_pds_a2_sample(i_duration, i_return_period, n_sample, flag)
% Return the sample of the rainfall distribution for a given 
%   duration and intensity based on Ann Arbor PDS data 
%   from NOAA Atlas 14. 
%
% The duration and return period are based on the PDS indices
%  eg. Return Per is: [0       1       2       5       10      25      50      100     200     500     1000] years
%        Duration is: [0 min   5min    10min   15min   30min   60min   2hr     3hr     6hr     12hr    24hr    2day...] 
%
%  A 10 year, 24 hour storm would be indexed (11,5)
%  --> The precip data is retrieved by calling rain_pds(11,5)
%
% n_sample corresponds to the index from the random distribution. 
%  To get multiple samples, insert an array (e.g. n_sample = 1:10)
%
% Flag options are: mean, lower, upper. These return the average rainfall
%   and the lower,upper values (resp.) for approx 1-sigma from the mean
%
% example calls:
% >>rain_pds_a2_sample(i_duration, i_return_period,19)
% ans =
%     3.7840
% 
% >>rain_pds_a2_sample(i_duration, i_return_period,[17 19])
% ans =
%     3.1963    3.7840 
%     
% >>rain_pds_a2_sample(i_duration, i_return_period,19,'upper')
% ans =
%     3.6500    

% Load Partial Duration Series (PDS) data for Ann Arbor
rain_pds_a2;

% Get the average and lower, upper 90-percentile.
%  Let the lower, upper percentiles approx the standard deviation
pds_avg  = rain_pds(i_duration, i_return_period);
pds_lower= rain_pds_lower(i_duration, i_return_period);
pds_upper= rain_pds_upper(i_duration, i_return_period);

% Represent the distribution as a log-normal distribution
%  Formulate the mean (m) and the variance (v) based on the pds data above
%  Back-calculate mu, sigma for the corresponding normal distribution
%    and use that as the input to lognstat()
% e.g. m = 1; v = 2;
m = pds_avg;
v = mean([ pds_lower - m
    m - pds_upper ])^2;

% Check the flag and adjust accordingly
if nargin > 3
    v = 0;
    if strcmpi(flag, 'mean')
        % pass
    elseif strcmpi(flag, 'lower')
        m = pds_lower;        
    elseif strcmpi(flag, 'upper')
        m = pds_upper;
    else
        warning('rain_pds_a2_sample.m: Invalid flag - %s. Returning mean rainfall..',flag)
    end
end

mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));

% Confirm the mean (M) and variance (V) of the resulting distribution
%  agrees with m, v from the PDS data
[M,V]= lognstat(mu,sigma);
% e.g.
% M =
%      1
% V =
%     2.0000

if abs( M-m ) > .01
    warning('rain_pds_a2_sample.m: Numerical imprecision -- abs( M-m ) > .01')
end

if abs( V-v ) > .01
    warning('rain_pds_a2_sample.m: Numerical imprecision -- abs( V-v ) > .01')
end


%%

% Seed the random number generator
rng(0)

% Generate the samples from the distribution and return the
%  result for the n_sample'th sample
n_rows = 1;
n_samples = max(n_sample);
X = lognrnd(mu,sigma,n_rows,n_samples);

inches_rainfall = X(n_sample);

% Generate a plot of the distribution
%{
X = lognrnd(mu,sigma,1,1e7);

MX = mean(X)
% MX =
%     0.9974
VX = var(X)
% VX =
%     1.9776

figure;hist(X,200)
xlim([0 12])
%}

% Generate 20 random samples from the distribution and
%  and verify the preceding samples are consistent
%{
for rnd_len = 1:20
    rng(0)
    
    X = lognrnd(mu,sigma,1,rnd_len);
    inches_rain = X(end);
    fprintf('%.4f: ', inches_rain); disp(X)
    
    clear X
end
%}