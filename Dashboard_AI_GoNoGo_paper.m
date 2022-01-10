%% Dashboard for simulating active inference for Go/No-go Avoid/Escape
% Exploring active-escape bias (Millner et al., 2019)

%% Single subject simulation

% first, set parameters for simulation
% -------------------------------------------------------------------------

% task parameters
mdp.n_trials = 200;   % the number of trials
mdp.p_Go     = 0.5;   % fraction of Go trials
mdp.p_av     = 0.5;   % fraction of avoid trials
mdp.rev_n    = 100;   % trial index when the reversal happens

% model parameters
mdp.c     = 8;        % sensitivity to the aversive outcome
mdp.k     = 0.1;      % stress weight for learning
mdp.alpha = 3;        % precision of action selection
mdp.beta  = 1;        % stochasticity of policy selection
mdp.e     = [1 1 1]'; % prior preference for policies: go, no-go, Pavlov
mdp.ups   = 20;       % the strength of the uniform prior for 'a' and 'b'
mdp.eta   = 1;        % base learning rate parameter
mdp.z     = 0.3;      % belief that Pavlovian policy leads to safety
mdp.w0    = 0.5;      % controllability threshold

% logistic decay parameter function parameters
mdp.df = 0;           % decay parameter (df=0 - SAPE-dependent, df>0 - fixed)
mdp.dm = 1.3;         % belief decay threshold
mdp.dl = 2;           % min
mdp.dh = 50;          % max
mdp.dg = 8;           % gradient

% run the simulation and plot results (Fig 4.)
% note, the results will vary across runs due to stochasticity
% -------------------------------------------------------------------------
MDP = GoNoGo_EA(mdp,'b');       % generate task structure
rng('shuffle');                 % pick a random seed 
MDP = spm_MDP_VB_LC_EA(MDP);    % simulate performance in the task

% set defaults and options for plotting
opts.mainplot = 1;
opts.showparams = {'k','dm','c','w0'};
set(0,'defaultAxesFontSize',12);    

plot_GoNoGo_EA_stats(MDP,opts); % plot results

% run the simulation and plot results (Fig 5.)
% note, the results will vary across runs due to stochasticity
% -------------------------------------------------------------------------
mdp.k = 1;                      % aversive learning scaling param
MDP = GoNoGo_EA(mdp,'b');       % generate task structure
rng('shuffle');                 % pick a random seed 
MDP = spm_MDP_VB_LC_EA(MDP);    % simulate performance in the task

% set defaults and options for plotting
set(0,'defaultAxesFontSize',12);    
opts.mainplot = 1;
opts.showparams = {'k','dm','c','w0'};

plot_GoNoGo_EA_stats(MDP,opts); % plot results

%% Simulate: vary one parameter at a time, but do so for many params
% (Fig. 6)

mdp.runs = 50;               % number of runs to average over

% set task parameters
mdp.n_trials = 400;          % the number of trials
mdp.rev_n    = 200;          % trial index when the reversal happens

mdp.learnwhat = 'b';         % which param to learn
opts.mainplot = 0;           % main plot for each run

opts.showparams = {'k','dm','c','w0','z'}; % display params on the plots

% set global param values
mdp.k  = 0.7;                % stress weight for learning
mdp.dm = 1.3;                % belief decay threshold 
mdp.c  = 8;                  % sensitivity to the aversive outcome
mdp.w0 = 0.6;                % controllability threshold
mdp.z  = 0.3;                % belief that Pavlovian policy leads to safety

% k
mdp.pn = {'k'};                     % choose params
mdp.pv = linspace(0.1,1.1,5);       % set a range of values
plot_GoNoGo_EA_many_many(mdp,opts)  % simulate and plot
    
% m
mdp.pn = {'dm'};                    % choose params
mdp.pv = linspace(0,2,5);           % set a range of values
plot_GoNoGo_EA_many_many(mdp,opts)  % simulate and plot

% w0
mdp.k = 0.9;
mdp.pn = {'w0'};                    % choose params
mdp.pv = linspace(0.2,0.8,5);       % set a range of values    
plot_GoNoGo_EA_many_many(mdp,opts)  % simulate and plot

% c
mdp.k  = 0.6;
mdp.dm = 1;
mdp.w0 = 0.5;
mdp.pn = {'c'};                     % choose params
mdp.pv = linspace(4,20,5);          % set a range of values
plot_GoNoGo_EA_many_many(mdp,opts)  % simulate and plot


%% Simulate: different subtypes
% (Fig. 7)

% set task parameters
mdp.n_trials = 400;          % the number of trials
mdp.rev_n    = 200;          % trial index when the reversal happens

% set global param values
mdp.k  = 0.6;                % stress weight for learning
mdp.dm = 1.3;                % belief decay threshold 
mdp.c  = 8;                  % sensitivity to the aversive outcome
mdp.w0 = 0.6;                % controllability threshold
mdp.z  = 0.3;                % belief that Pavlovian policy leads to safety

% set defaults and options for plotting
opts.showparams = {'k','dm','c','w0'};
set(0,'defaultAxesFontSize',12); 
opts.mainplot = 0;

% subtype m
mdp.k  = 0.7;                    % stress weight for learning
mdp.dm = 2;                      % belief decay threshold 
mdp.c  = 8;                      % sensitivity to the aversive outcome
mdp.w0 = 0.6;                    % controllability threshold

MDP = GoNoGo_EA(mdp,'b');        % generate task structure
rng('shuffle');                  % pick a random seed 
MDP = spm_MDP_VB_LC_EA(MDP);     % simulate performance in the task
plot_GoNoGo_EA_stats(MDP,opts);  % plot results

% subtype w0
mdp.k  = 0.9;                    % stress weight for learning
mdp.dm = 1.3;                    % belief decay threshold 
mdp.c  = 8;                      % sensitivity to the aversive outcome
mdp.w0 = 0.8;                    % controllability threshold

MDP = GoNoGo_EA(mdp,'b');        % generate task structure
rng('shuffle');                  % pick a random seed 
MDP = spm_MDP_VB_LC_EA(MDP);     % simulate performance in the task
plot_GoNoGo_EA_stats(MDP,opts);  % plot results

% subtype k
mdp.k  = 1.1;                    % stress weight for learning
mdp.dm = 1.3;                    % blief decay threshold 
mdp.c  = 8;                      % sensitivity to the aversive outcome
mdp.w0 = 0.6;                    % controllability threshold

MDP = GoNoGo_EA(mdp,'b');        % generate task structure
rng('shuffle');                  % pick a random seed 
MDP = spm_MDP_VB_LC_EA(MDP);     % simulate performance in the task
plot_GoNoGo_EA_stats(MDP,opts);  % plot results

% subtype c
mdp.k  = 0.6;                    % aversive learning scaling param
mdp.dm = 1;                      % belief decay threshold 
mdp.c  = 20;                     % sensitivity to the aversive outcome
mdp.w0 = 0.5;                    % controllability threshold

MDP = GoNoGo_EA(mdp,'b');        % generate task structure
rng('shuffle');                  % pick a random seed 
MDP = spm_MDP_VB_LC_EA(MDP);     % simulate performance in the task
plot_GoNoGo_EA_stats(MDP,opts);  % plot results
