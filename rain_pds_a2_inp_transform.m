function inp = rain_pds_a2_inp_transform(n_sample, i_duration, i_return_period, flag)
% Requires
%  template = [pwd '/runfile_xw_a2_template.inp'];
%
% Set the duration and return period based on the PDS indices
%  eg. Return Per is: [0       1       2       5       10      25      50      100     200     500     1000] years
%        Duration is: [0 min   5min    10min   15min   30min   60min   2hr     3hr     6hr     12hr    24hr    2day...] 
%
%  A 10 year, 24 hour storm would be indexed (11,5)
%  --> The precip data is retrieved by calling rain_pds(11,5)
%
if nargin < 2
    i_duration      = 11; % 11 corresponds to 24 hr storm
    i_return_period = 5;  % 5  corresponds to 10 yr return period
end

% Load Partial Duration Series (PDS) data for Ann Arbor
rain_pds_a2;

%{
fprintf('%15s: %ghr\n%15s: %gyr\n',...              % Print the storm info
    'Duration',      rain_pds(i_duration,1)/60, ... % Print the duration rain_pds(_,1)
    'Return Period', rain_pds(1,i_return_period))   % Print the return period rain_pds(1,_)
%}

% Define timestep, duration, and total rainfall for generating the design storm
dt = 5/60;                      % 5 mins in hours
dur= rain_pds(i_duration,1)/60; % hours

if nargin == 4
    total_in_rainfall = rain_pds_a2_sample(i_duration, i_return_period, n_sample, flag);
else
    total_in_rainfall = rain_pds_a2_sample(i_duration, i_return_period, n_sample);
end
% total_in_rainfall = rain_pds(i_duration, i_return_period); % need to first call rain_pds_a2 


% Get the timeseries for an SCS II design storm
%  This returns the cumulative rainfall
%  Max recommended storm length is 24 hours
[rain_transform,t]=scs_II_transform(dt,dur,total_in_rainfall);

% Convert cumulative rainfall to intensity by taking the difference 
%  between each datapoint
rain_transform = [0 diff(rain_transform) 0];
t = [t t(end) + t(end)-t(end-1)];

% Separate the hours and mins from the time array orginally in decimal hours
hours = floor(t);
t_mins = t - hours;
mins   = round(t_mins * 60);


% Append timeseries to prepared .inp-file template
%  Define the filepaths
%inp = [pwd sprintf('/runfile_xw_a2_template_%06g.inp',n_sample)];
template = [pwd '/runfile_xw_2_a2_template.inp'];
inp = [pwd sprintf('/runfile.inp')];

%  Create the file
copyfile(template, inp);
fid = fopen(inp,'a');

%fwrite(fid, sprintf('\n'));

% Append timeseries to prepared .inp-file template
for m = 1:length(rain_transform)
    %fprintf('%s	          	%02d:%02d:%02d  	%f\n','DESIGN_10YR12HR_ALT',hours(m),mins(m),0,rain_transform(m))
    tmp_line = sprintf('%s	          	%02d:%02d:%02d  	%f\n','DESIGN_10YR12HR_ALT',hours(m),mins(m),0,rain_transform(m));
    fwrite(fid,tmp_line);
end

fclose(fid);


if abs( total_in_rainfall - trapz(t/dt,rain_transform) ) > .01
    warning('rain_pds_a2_inp_transform.m: Numerical imprecision for rain_pds(%g,%g)\n',i_duration,i_return_period)
end
cumtrapz(t / dt,rain_transform); % check cumulative rainfall