function [] = swmm_tracking(i_combination, n_sample, input_file)

pct_full = .9;
pct_full_final = .9;

n_controllers = 11;
plot_response = 1;
save_plot     = 0;

%% File Setup
%

if nargin < 3
    % Requires "runfile.inp" exists in current directory
    % even though runfile_xw_2_a2_template.ing is used instead
    %rain_pds_a2_inp_transform(n_sample);
    
    inp = 'runfile.inp';
    inp_name = '';
else
    warning('swmm_tracking.m: Not Tested. Exiting...')
    return 
    
    inp = input_file;
    [~,inp_name,~] = fileparts(inp);
    inp_name = [inp_name '/'];
end

% copy the appropriate library as swmm.so depending on the operating system
if ispc
    copyfile('swmm5.dll','swmm5.so')
elseif ismac
    copyfile('swmm5_osx.so','swmm5.so')
else % isunix
    copyfile('swmm5_linux.so','swmm5.so')
end

sprintf('sim/%s%06g/%06g/',inp_name,i_combination,n_sample);
files = {'create_system_matrices_v2.m', 'SWMM.m', 'read_orifice_from_node.m', 'swmm2lqr_minreal_v4.m',...
    'read_storage_curves.m', 'swmm5.h', 'create_labels.m', 'read_subcatch_outlet.m', 'swmm5.so', ...
    'create_storage_block_v2.m', 'swmm_properties.mat', 'plot_swmm2lqr.m', 'reachable_test.m',...
    'rain_pds_a2_inp_transform.m','rain_pds_a2.m','rain_pds_a2_sample.m','scs_II_transform.m',...
    'runfile_xw_2_a2_template.inp'};
files{end+1} = inp;

% make the directories to save the data and errors
OUTDIR  = sprintf('data/%s%06g/',inp_name, i_combination);
mkdir(OUTDIR);

% make the directory for the simulation files
DIR = sprintf('sim/%s%06g/%06g/', inp_name, i_combination, n_sample);
mkdir( DIR )

% copy all files to i_combination
for i_file = 1:length(files)
    copyfile(files{i_file}, [DIR files{i_file}])
end

% change directory to i_combination
cd(DIR)

% Check if simulation has already been completed
if nargin < 3
    OUTFILE = sprintf('../../../data/%06g/%06g.mat',i_combination,n_sample);
    rain_pds_a2_inp_transform(n_sample,11,5,'mean');
    inp = 'runfile.inp';
else
    OUTFILE = sprintf('../../../../data/%s%06g/%06g.mat',inp_name,i_combination,n_sample);
end
if exist(OUTFILE,'file') ~= 0
% (change to n_sample if the combination stays the same but you'd like to
%  vary the rainfall intensity)    
    fprintf('File %s already exists! Simulation %g, %g skipped.\n',OUTFILE,i_combination,n_sample);
    return
end

%inp = 'runfile.inp';

% Load the properties of the model -- this must be saved a priori for HPC
load('swmm_properties')

%% Load the SWMM library
%

swmm = SWMM;


%%
% Define function handles

clip = @(x, x_min, x_max) max( min( x, x_max), x_min);
remove_blanks = @(array) array(~ismember(array,{'','(null)'}));

%% LQR Setup

metric  = 0;
ft2m    = 0.3048;

rho     = 35; % Ratio of cost on states vs inputs
q       = 100 * ones(1,11); % Cost on states
Q_min   = zeros(1,11) * 0.3048^3;
Q_max   = 10 * ones(1,11)  * 0.3048^3; % cfs to m3/s

Q_max   = Q_max * 2; % Apparently I need a factor of 2 somewhere

%% Set up the controllable vs uncontrollable orifices
%

% orifice properties
g = 9.81;   % gravity [m/s]
Cg = 1.00;   % gate coeff calibration [-]
ug = 0.63;   % gate coeff [-]
Wg = 0.3048; % gate width [-]
h2 = 0;      % setpoint for submerged height downstream [m]
hcr= 0;      % height of crest underneath gate [m]

% estimate gate opening by inverting flow through a submerged gate
A_gate = @(Qc,h1)    real( 1./(Cg*ug * sqrt(2*g*(h1 - h2))) ) .* Qc ;
lambda = @(Qc,h1,Wg) A_gate(Qc,h1) ./ (Wg .* Qc);
A_gate_clipped = @(Qc,h1,Amin,Amax)    clip(A_gate(Qc,h1), Amin, Amax);
lambda_clipped = @(Qc,h1,Amin,Amax,Wg) A_gate_clipped(Qc,h1, Amin, Amax) ./ (Wg .* Qc);

% 0.3605 <-- g*Cg*Wg*ug*hg(1)/sqrt(2*g*(h1(1)-ug*hg(1))) where Cg = 2.5
% length(storages) == 11
lambda_array = 0.3605 * ones(1,11);
G_lambda_array = 0.3605 * ones(1,11);

% Construct an array of all possible controller combinations
n_orifices = length(lambda_array);
A_ctl = cell(n_controllers,1);
v = 1:n_orifices;

for n = 1:n_controllers
    A_ctl{n} = nchoosek(v,n);
end

ncol   = cellfun('size',A_ctl,2);
maxcol = n_orifices;
for n = 1:length(A_ctl)
    if ncol(n) < maxcol
        A_ctl{n}(end,maxcol) = 0;
    end
end

combinations = cat(1, A_ctl{:});

% Set the controlled valves to NaN
if ismember(i_combination, 1:length(combinations))
    % Check if i_combination is a valid permutation index. Otherwise, default to no control
    lambda_array( combinations( i_combination , combinations(i_combination,:)>0) ) = nan;
else
    fprintf('Default no control: i_combination = %d and there are 1:%d possible indices\n',i_combination, length(combinations));
end

%% SWMM simulation parameters
%

% Load the properties of the model
%{
[I,nodes,links]=swmm.get_incidence_matrix(inp);

storages = remove_blanks( swmm.get_all(inp,swmm.STORAGE,swmm.NONE) );
junctions= remove_blanks( swmm.get_all(inp,swmm.JUNCTION,swmm.NONE) );
conduits = remove_blanks( swmm.get_all(inp,swmm.LINK,swmm.NONE) );
orifices = remove_blanks( swmm.get_all(inp,swmm.ORIFICE,swmm.NONE) );
subcatch = remove_blanks( swmm.get_all(inp,swmm.SUBCATCH,swmm.NONE) );
outfalls = nodes( ~ismember(nodes,union(storages,junctions)) );


if length(lambda_array)   ~= length(orifices)
    error('lambda_array: Incorrect number of passive inputs');
end

if length(G_lambda_array) ~= length(orifices)
    error('G_lambda_array: Incorrect number of control inputs');
end

if length(Q_max) ~= length(orifices) || length(Q_min) ~= length(orifices)
    error('Q_max, Q_min: Incorrect number of min/max pairs')
end

if length(q) ~= length(storages)
    error('q: Incorrect number of weights for pond heights')
end
%}

%%% Use the above properties to
% - estimate storage curves,
% - associate outlets with subcatchments
read_storage_curves
read_subcatch_outlet
read_orifice_from_node

if ~metric
    ft2m = 0.3048; % 1 ft = 0.3048 m
    
    % Assume all three Maps have the same key
    keys = storage_curves_A.keys;
    
    for k = 1:length(keys)
        
        storage_curves_h(keys{k})  = storage_curves_h(keys{k}) * ft2m;
        storage_curves_A(keys{k})  = storage_curves_A(keys{k}) * ft2m^2;
        storage_curves_As(keys{k}) = storage_curves_As(keys{k}) * ft2m^2;
        
    end
end

% figure_path = sprintf('%s%s%02d controllers/storm %05d/',pwd,'/swmm_files/',1, nstorm );

%%% Run the simulation without any control to get the rainfall, runoff,
% as well as estimate peak depth and flow

warning off

% BIG NOTE: subcatch, storages, orifices do NOT necessarily return in
% corresponding order (ie depth(:,i) may really correspond with flow(:,j)
[sim_error, sim_duration] = swmm.run_simulation(inp);
[t, runoff] = swmm.read_results(subcatch,swmm.SUBCATCH,swmm.RUNOFF);
[~, precip] = swmm.read_results(subcatch,swmm.SUBCATCH,swmm.PRECIPITATION);
[~, depth ] = swmm.read_results(storages,swmm.NODE,swmm.DEPTH); % height, SI
[~, flow  ] = swmm.read_results(orifices,swmm.LINK,swmm.FLOW); % flow, SI

warning on

%t = [0:1/12:126]'; % Assume 126 hr runtime at 5/60 = 1/12 hr intervals
%T = 300; %5 minutes * 60s/min -- the timestep between each control input

t_seconds   = t * 3600;
T           = (t(2)-t(1)) * 3600;

%%
% Added in case we have to use a configuration w/o control
% - Save the variables from the initial swmm model run and return
% - This code is copy/pasted from below

if ~ismember(i_combination, 1:length(combinations))
    [t, tmp_runoff] = swmm.read_results(subcatch,swmm.SUBCATCH,swmm.RUNOFF);
    [~, tmp_precip] = swmm.read_results(subcatch,swmm.SUBCATCH,swmm.PRECIPITATION);
    [~, tmp_flood ] = swmm.read_results(storages,swmm.NODE,swmm.FLOODING); % height, SI
    [~, tmp_depth ] = swmm.read_results(storages,swmm.NODE,swmm.DEPTH); % height, SI
    [~, tmp_flow  ] = swmm.read_results(orifices,swmm.LINK,swmm.FLOW); % flow, SI
    [~, flood ] = swmm.read_results(storages,swmm.NODE,swmm.FLOODING);
    
    try
        n_tsteps = ceil( t(end)/dt );
        OR_setting   = ones( n_tsteps, length(G_lambda_array) );
        valve  = OR_setting;
    catch
    end
    
    save(OUTFILE)
    return
end


%% Initialize additional hash maps used for system dynamics
%

% Create a hash-map of the lambda values with the key as the name of the
% corresponding storage node
lambdas   = containers.Map(storages,lambda_array);   % valves that are 100% open
G_lambdas = containers.Map(storages,G_lambda_array); % valves that are controllable by matrix G
heights   = containers.Map(storages, pct_full_final*3.048+zeros(size(lambda_array)) ); % storage heights, initialized at zero meters

%% Linearized Dynamics Observations and Control Matrices
% Create matrices for xdot = Ax+Bd+Gu; y = Cx + Du
% - A is the dynamics matrix
% - B is the disturbance matrix (for runoff)
% - C is the observation matrix
% - D is the feedthrough matrix
% - G is the control matrix

show_component = 0;
create_system_matrices_v2

%% LQR Controller Design
% Calculate the gain, K, and close the loop given cost matrices Q, R
% - Q is the cost on state variables
% - R is the cost on input actions
%
% Now we write the system dynamics as xdot = (A+GK)x + Bd

% Isolate the rows that correspond to pond heights
Q_cols = C_h;%G > 0;
Q_cols = Q_cols * diag( sqrt(q) );
Q = rho * (Q_cols*(Q_cols'));

r = ones(size(lambda_array));
R = diag( r(isnan(lambda_array)) );

% Calculate the gain, K
%
% Returns: U, Ahat, Bhat, Chat, Dhat, Qhat, Rhat, Khat
% where Bhat is from the gain matrix, -G (NOTE: T/As was multiplied into G in v2)
%
% Note:
% Atil =  U*A*U'; %Atil = U*A*inv(U); Atil = U*A/U;
% Btil = -U*G;
% Ktil =  K*U';
swmm2lqr_minreal_v4;

% Posterity -- Calculate the gain, K
%   K = dlqr(A,-T/As(1)*G, Q, R); % Posterity
% Posterity -- Create the closed loop system
%   sys_ctl = ss(A+B*K,B,C,D,T);
%   sys_ctl = ss(A+T/As(1)*G*K, B,C,D,T);


%% Simulate the open-loop model and compare to SWMM
%
%{
% create open-loop model (no control)
sys = ss(A,B,C,D,T);

n_storage = length(storages);

% Create input vectors for lsim()
u = zeros(length(t), size(B,2));
%if length(subcatch) < n_storage
    u(:,1:length(subcatch)) = runoff;
%end
u      = u * 0.0283168; % convert Q_runoff: cfs to cubic meter/s

%%% Simulate the open-loop dynamics %%%

    [y_metric,~,x] = lsim(sys,u,t_seconds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert from metric to US
% * 0.3048 ft       = 1 meter
% * 0.0283168 cu ft = 1 cu. m
m2ft = ones(1,n_storage) * 1/0.3048;
y = y_metric * diag([m2ft m2ft.^3]);

OR_setting = 0*y(:,1:n_storage)+1;
%}
% Plot response
%{
if plot_response
    f = figure('units','normalized','outerposition',[0 0 1 1]); %figure;
    
    % Get the indices to associate the correct orifices with their storage
    % nodes. The call to Matswmm above generates the storages and orifices
    % in ABC order, which isn't always the case
    [~,i_sort] = sort(orifice_from_node.values);
    [~,i_rev]  = sort(i_sort);
    
    ax(1) = subplot(4,1,1);
    bar(t,u);
    set(ax(1),'Ydir','reverse')
    
    ax(2) = subplot(4,1,2);
    plot(t,y(:,1:n_storage),'-x','linewidth',2); %ylim([-1 7.5])
    hold on;plot(t,depth,'linewidth',4);%ylim([-1 7.5])
    
    height_labels = create_labels(storages,'h');
    legend(height_labels)
    
    ax(3) = subplot(4,1,3);
    plot(t,y(:, 1+n_storage:end),'-x','linewidth',2);
    hold on;plot(t,flow,'linewidth',4);%ylim([-1 7.5])
    
    flow_labels = create_labels(storages,'Q');
    legend(flow_labels)
    
    ax(4) = subplot(4,1,4);
    plot(t,OR_setting,'linewidth',2);
    legend(orifices(i_rev))
    
    linkaxes(ax,'x')
    xlim([0 10*round(max(t)/10)+10])
    
    mkdir(figure_path)
    %hgexport(f,[figure_path 'open_loop'])
    export_fig([figure_path 'open_loop.png'],'-nocrop','-transparent','-r150',f)
end
%}
%% Simulate the closed-loop model with clipping and compare to SWMM
%
%{
n_storage = length(storages);

m2ft = ones(1,n_storage) * 1/0.3048;

% Create input vectors for lsim()
u = zeros(length(t), size(B,2));
%if length(subcatch) < n_storage
    u(:,1:length(subcatch)) = runoff;
%end
u      = u * 0.0283168; % convert Q_runoff: cfs to cubic meter/s

%%% Simulate the open-loop dynamics %%%

    %[y_metric,~,x] = lsim(sys,u,t_seconds);
    
    xpp = zeros( length(t), size(A,2) );
    y_metric = zeros( length(t), size(C,1) );
    
    Qout         = 0*r(isnan(lambda_array))';
    Qout_clipped = 0*r(isnan(lambda_array))';
    
    Qout_min     = Q_min(isnan(lambda_array))';
    Qout_max     = Q_max(isnan(lambda_array))';
    
    max_lambda = G_lambda_array(isnan(lambda_array));
    
    xd = sum(C_h,2)*[9+11/12]*ft2m;
    
    for step = 1:length(t)
       
        Qout(:,step) = (-K*(xpp(step,:)' - xd))';    % metric, [m^3 s^-1]
        
        % Try the minimal realization controller
        xtil = U*(xpp(step,:)' - xd - 0); %replace "- 0" with mean rain forecast
        xhat = xtil(1:size(Ahat));
        Qout(:,step) = (-Khat*xhat)';
        
        Qout_clipped(:,step) = clip( Qout(:,step), Qout_min, Qout_max ); % metric, [m^3 s^-1]
        
        heights   = containers.Map(storages, y_metric(step,1:length(storages) ));
        y(step,:) = y_metric(step,:) * diag([m2ft m2ft.^3]);
        
        % Compute area to open gate based on desired outflow and current height
        
        tmp = (Qout_clipped(:,step) * 1/0.0283168)'./y(step,find([isnan(lambda_array) 0*isnan(lambda_array)]))* 0.3048^2;
        %tmp = (Qout(:,step) * 1/0.0283168)'./y(step,[1 3])* 0.3048^2;
        tmp( tmp < 0 | isnan(tmp) ) = 0;
        tmp( tmp > max_lambda ) = max_lambda( tmp > max_lambda);
        
        lambda_ctl(step,:) = tmp;
        
        tmp = lambda_array;
        tmp(isnan(lambda_array)) = lambda_ctl(step,:);
        lambdas = containers.Map(storages,tmp);
        
        
        % Save the orifice settings for linearized case
        OR_setting(step,:) = tmp ./ G_lambda_array;
        
        create_system_matrices_v2
        % Try just updating G and using the second xpp equation
        
        xpp(step+1,:) = A*xpp(step,:)' + B*u(step,:)';
        %xpp(step+1,:) = (A+T/As(1)*G*K)*xpp(step,:)' + B*u(step,:)';
        y_metric(step+1,:) = C*xpp(step+1,:)';
        
    end

x = xpp(2:end,:);
y_metric = y_metric(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert from metric to US
% * 0.3048 meter    = 1 ft
% * 0.0283168 cu. m = 1 cu ft
m2ft = ones(1,n_storage) * 1/0.3048;
y = y_metric * diag([m2ft m2ft.^3]);
%}
% Plot response
%{
if plot_response
    f = figure('units','normalized','outerposition',[0 0 1 1]); %figure;
    
    ax(1) = subplot(4,1,1);
    bar(t,u(:,1));
    set(ax(1),'Ydir','reverse')
    
    ax(2) = subplot(4,1,2);
    plot(t,y(:,1:n_storage),'-x','linewidth',2); %ylim([-1 7.5])
    %hold on;plot(t,depth,'linewidth',4);%ylim([-1 7.5])
    
    height_labels = create_labels(storages,'h');
    legend(height_labels)
    
    ax(3) = subplot(4,1,3);
    plot(t,y(:, 1+n_storage:end),'-x','linewidth',2);
    %hold on;plot(t,flow,'linewidth',4);%ylim([-1 7.5])
    
    flow_labels = create_labels(storages,'Q');
    legend(flow_labels)
    
    ax(4) = subplot(4,1,4);
    plot(t,OR_setting,'linewidth',2);
    legend(orifices(i_rev))
    
    linkaxes(ax,'x')
    xlim([0 10*round(max(t)/10)+10])
    
    OR_setting_linear = OR_setting;
    y_linear = y;
    
    %hgexport(f,[figure_path 'closed_loop_linear'])
    export_fig([figure_path 'closed_loop_linear.png'],'-nocrop','-transparent','-r150',f)
end
%}
%%
clear OR_setting y

t_start = tic;
m2ft = ones(1,n_storage) * 1/0.3048;

swmm.initialize(inp);

%%% Simulate the closed-loop dynamics %%%

dt = 1/120; % Assume MIMINUM_STEP = 0.5 (or 30 seconds)
n_lambda = sum( isnan(lambda_array) );
n_tsteps = ceil( t(end)/dt );

sparse_A = sparse(A);
sparse_B = sparse(B);

xpp = zeros(  n_tsteps, size(A,2) );
y_metric = zeros( n_tsteps, size(C,1) );
y        = zeros( n_tsteps, size(C,1) );
u_tmp    = zeros( 1,size(B,2) );

Qout         = zeros( n_lambda, n_tsteps ); % 0*r(isnan(lambda_array))';
Qout_clipped = zeros( n_lambda, n_tsteps ); % 0*r(isnan(lambda_array))';

lambda_ctl   = zeros( n_tsteps, n_lambda );
OR_setting   = zeros( n_tsteps, length(G_lambda_array) );

depth_swmm   = zeros( n_tsteps, size(C_h',1) ); % Recall C_* is transposed, use size(*,1) for consistency
flow_swmm    = zeros( n_tsteps, size(C_Q',1) );
t_swmm       = zeros( n_tsteps, 1);

Qout_min     = Q_min(isnan(lambda_array))';
Qout_max     = Q_max(isnan(lambda_array))';

max_lambda = G_lambda_array(isnan(lambda_array));

xd = 0;
xref = pct_full_final*3.048 + zeros(n_mv,1);
%xref = pct_full_final*3.048 * sum(C_h>0,2);

% Get the indices to associate the correct orifices with their storage
% nodes. The call to Matswmm above generates the storages and orifices
% in ABC order, which isn't always the case
[~,i_sort] = sort(orifice_from_node.values);
[~,i_rev]  = sort(i_sort);

%     OR_setting = 0*q(:);
%     lambda_ctl = 0*q(:);

step = 1;

%control_sequence = reachable_test(A,B,C_h,D,G,runoff,pct_full,pct_full_final);
% Use the control sequence from runfile_xw_a2_10yr_24hr
control_sequence = ...
    [-0.0128    0.0155    0.0249    0.0411    0.0114    0.0042   -0.0061    0.0908    0.0875    0.0010    0.0995];
    % worse --> [0         0    0.1278         0    0.0329    0.0115         0    0.4464    0.5621         0    0.4871];

while ~swmm.is_over
    
    % Set the gate position using gain matrix from LQR
    
    %Qout(:,step) = (-K*(xpp(step,:)' - xd - 0))';    % metric, [m^3 s^-1]
    
    % Try the minimal realization controller
    xtil = U*xpp(step,:)';     % metric, [m^3 s^-1]
    xhat = xtil(1:size(Ahat));
    Qout(:,step) = (-Khat*xhat+E*xref)'; %(-Khat*xhat);
    %Qout(:,step) = (-Khat*(xhat-xref))'; %(-Khat*xhat);
    
    
    % Clip to limit the outflow
    Qout_clipped(:,step) = clip( Qout(:,step), Qout_min, Qout_max ); % metric, [m^3 s^-1]
    
    y(step,:) = y_metric(step,:) * diag([m2ft m2ft.^3]);
    
    %{
    % update the dynamics and gain matrices every 5 mins (or 60 steps)
    if ~mod(step,60)
        
        %fprintf('%g\n',t_swmm(step-1))
        
        heights   = containers.Map(storages, y_metric(step,1:size(y_metric,2)/2) ); % storage heights, initialized at zero meters
        create_system_matrices_v2
        swmm2lqr_minreal_v4
        
        %pct_full = y(step,1:size(y_metric,2)/2) / 10; %pct for 10ft full
        %control_sequence = reachable_test(A,B,C_h,D,G,runoff,pct_full,pct_full_final);
        
    end
    %}
    
    % Convert flow to orifice area; A_gate is constrained below by 0 and above by max gate area
    A_min = 0;
    A_max = 2*Wg*Wg;
    A_gate_k = A_gate_clipped( Qout_clipped(:,step)', y_metric(step,find([isnan(lambda_array) 0*isnan(lambda_array)])) , A_min, A_max);
    
    % Convert Area to percentage opening to plug into MatSwmm
    OR_setting(step, :)                   = ones( size(G_lambda_array) );
    OR_setting(step, isnan(lambda_array)) = A_gate_k ./ A_max;
    
    % Calculate lambda corresponding to the clipped gate area
    % Assume all gates at the same width, Wg
    % Recall, h_gate = lambda * Q_clipped
    % --> A_gate = h_gate * Wg
    % --> A_gate = (lambda * Q_clipped) * Wg
    % --> A_gate = lambda * (Q_clipped * Wg)
    % --> lambda = A_gate / (Q_clipped * Wg)
    
    lambda_ctl(step,:) = lambda_clipped( Qout_clipped(:,step)', y_metric(step,find([isnan(lambda_array) 0*isnan(lambda_array)])) , A_min, A_max, Wg);
    
    %%% Old conversion from flow to orifice area
    % tmp = (Qout_clipped(:,step) * 1/0.0283168)'./y(step,find([isnan(lambda_array) 0*isnan(lambda_array)]))* 0.3048^2;
    % %tmp = (Qout(:,step) * 1/0.0283168)'./y(step,[1 3])* 0.3048^2;
    % tmp( tmp < 0 | isnan(tmp) ) = 0;
    % tmp( tmp >  max_lambda) = max_lambda( tmp > max_lambda );
    %
    % lambda_ctl(step,:) = tmp;
    % OR_setting(step,:) = G_lambda_array;
    % OR_setting(step, isnan(lambda_array) ) = tmp;
    % OR_setting(step,:)                     = OR_setting(step,:) ./ G_lambda_array;
    % % Use G_lambda_array as the max lambda value to normalize the gates
    
    
    % Reorder the OR_setting from the order of the storage nodes
    % to now correspond with the order of the orifices and then
    % modify the settings in SWMM
    %swmm.modify_settings( orifices(i_rev), OR_setting(step,:) );
    swmm.modify_settings( orifices, OR_setting(step,:) );
    %swmm.modify_settings( orifices(isnan(lambda_array)), OR_setting(step,isnan(lambda_array)) );
    
    
    % Step forward
    t_swmm(step)  = swmm.run_step;
    
    depth_swmm(step,:) = swmm.step_get_fast(storages,swmm.DEPTH,swmm.SI);
    flow_swmm(step,:)  = swmm.step_get_fast(orifice_from_node.values,swmm.FLOW,swmm.SI);
    %flow_swmm(step,:)  = swmm.step_get(orifices,swmm.FLOW,swmm.SI);
    
    % step flows forward and add in runoff from subcatchment(s)
    u_tmp = u_tmp * 0;
    u_tmp(1:length(subcatch)) = swmm.step_get_fast(subcatch,swmm.RUNOFF,swmm.SI);
    
    %         if length(subcatch) > n_storage
    %             % NOTE: n_storage == size(B,2)
    %             tmp = tmp(1:n_storage);
    %         end
    
    % Assign flow & propagate down
    %xpp(step+1,:) = sparse_A*xpp(step,:)' + sparse_B*u_tmp';
    xpp(step+1,:) = 0*(A*xpp(step,:)' + B*u_tmp');
    
    % assign latest depths
    %  sum(C_h>0,2) returns a vector where a 1 corresponds to a depth
    %  state
    xpp(step+1, find(sum(C_h>0,2)) ) = depth_swmm(step,:);
    
    
    % update observations
    %         y_metric(i+1,:) = [
    %             swmm.get('S1',swmm.DEPTH,swmm.SI) % height
    %             xpp(i+1,n_delay+1)                % inflow
    %             swmm.get('OR1',swmm.FLOW,swmm.SI) % outflow
    %             ];
    
    obs_tmp = [depth_swmm(step,:) flow_swmm(step,:)]';
    y_metric(step+1,:) = obs_tmp(:);
    
    step = step+1;
end

t_swmm(end) = 2*t_swmm(end-1)-t_swmm(end-2); % = time(end-1) + diff( time(end-1), time(end-2) )



x = xpp(2:end,:);
y_metric = y_metric(2:end,:);

[errors, d] = swmm.finish();
t_end=toc;

[t, tmp_runoff] = swmm.read_results(subcatch,swmm.SUBCATCH,swmm.RUNOFF);
[~, tmp_precip] = swmm.read_results(subcatch,swmm.SUBCATCH,swmm.PRECIPITATION);
[~, tmp_flood ] = swmm.read_results(storages,swmm.NODE,swmm.FLOODING); % height, SI
[~, tmp_depth ] = swmm.read_results(storages,swmm.NODE,swmm.DEPTH); % height, SI
[~, tmp_flow  ] = swmm.read_results(orifices,swmm.LINK,swmm.FLOW); % flow, SI
[~, flood ] = swmm.read_results(storages,swmm.NODE,swmm.FLOODING);

try
    valve  = OR_setting;
catch
end

save(OUTFILE)

%%
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert from metric to US
% * 0.3048 meter    = 1 ft
% * 0.0283168 cu. m = 1 cu ft
m2ft = ones(1,n_storage) * 1/0.3048;
y = y_metric * diag([m2ft m2ft.^3]);
%}
% Plot response
if plot_response
    f = figure('units','normalized','outerposition',[0 0 1 1]); %figure;
    
    ax(1) = subplot(4,1,1);
    bar(t,runoff(:,1));
    set(ax(1),'Ydir','reverse')
    
    ax(2) = subplot(4,1,2);
    plot(t_swmm(1:50:end),y(1:50:end,1:n_storage),'-x','linewidth',2); %ylim([-1 7.5])
    %hold on;plot(t,depth,'linewidth',4);%ylim([-1 7.5])
    
    height_labels = create_labels(storages,'h');
    legend(height_labels)
    
    ax(3) = subplot(4,1,3);
    plot(t_swmm(1:50:end),y(1:50:end, 1+n_storage:end),'-x','linewidth',2);
    %hold on;plot(t,flow,'linewidth',4);%ylim([-1 7.5])
    
    flow_labels = create_labels(storages,'Q');
    legend(flow_labels)
    
    ax(4) = subplot(4,1,4);
    plot(t_swmm(1:50:end),OR_setting(1:50:end,:),'linewidth',2);
    legend(orifices(i_rev))
    
    linkaxes(ax,'x')
    xlim([0 10*round(max(t)/10)+10])
    
    if save_plot
        %hgexport(f,[figure_path 'closed_loop_swmm'])
        export_fig([figure_path 'closed_loop_swmm.png'],'-nocrop','-transparent','-r150',f)
    end
end
%%
%

%{
if plot_response
    plot_swmm2lqr
end
%}