function MDP = GoNoGo_EA(mdp,learnwhat)

rng('default'); % ensure the same trial sequence for all runs

%-------------------------------------------------------------------------
% This code sets task structure for active inference model for implementing
% Millner et al., (2019) Go/No-Go task with escape and avoid conditions 
% with an aversive stimulus.
%
% The agent starts at a neutral location 1 and makes a compulsory 'move' 
% to location 2, where either only a cue is presented (avoid condition) 
% or an outcome and a cue (escape condition). Go/No-Go responses take the
% agent to location 3 and 1 respectively. Correct responses result in a 
% neutral outcome in avoid condition and a discontinuation of the aversive 
% stimulus in escape condition. 

% Inputs: params   - structure containing parameter values
%         learwhat - assume learning in 'a' or 'b'
%
% Output: a completed MDP structure

% States:
%   1 = at location 1, GA context
%   2 = at location 1, NGA context
%   3 = at location 1, GE context
%   4 = at location 1  NGE context
%   5 = at location 2, GA cue
%   6 = at location 2, NGA cue
%   7 = at location 2, GE cue
%   8 = at location 2  NGE cue
%   9 = at location 3, sound off
%   10 = at location 3, sound on
%   11 = at location 1, sound off
%   12 = at location 1, sound on

% Observations:
%   1 = at 1 (no info on context)
%   2 = at 2, cue = GA
%   3 = at 2, cue = NGA
%   4 = at 2, cue = GE
%   5 = at 2, cue = NGE
%   6 = at 3, sound off
%   7 = at 3, sound on
%   8 = at 1, sound off
%   9 = at 1, sound on
%--------------------------------------------------------------------------

% Unpack relevant input parameters
n_trials = mdp.n_trials; 
p_Go     = mdp.p_Go;
p_av     = mdp.p_av; 
rev_n    = mdp.rev_n;
c        = mdp.c;
k        = mdp.k;
df       = mdp.df;
ups      = mdp.ups;
e        = mdp.e; 
z        = mdp.z;

% A - outcome probabilities given the state
%--------------------------------------------------------------------------
% Most outcomes are deterministic, apart from 20% chance of aversive
% outcome even with the correct go/no-go response if learning in A is
% assumed
%--------------------------------------------------------------------------

% If learning in A is assumed or not
if strcmp(learnwhat,'a')
    a1 = 1; a2 = 1 - a1;     % process likelihoods 
    a3 = 0.8; a4 = 1 - a3;   % initial model likelihoods
else
    a1 = 1; a2 = 1 - a1;     % process likelihoods 
    a3 = 1; a4 = 1 - a3;     % initial model likelihoods
end

A_ENV = [1	1	1	1	0	0	0	0	0	0	0	0;
         0	0	0	0	1	0	0	0	0	0	0	0;
         0	0	0	0	0	1	0	0	0	0	0	0;
         0	0	0	0	0	0	1	0	0	0	0	0;
         0	0	0	0	0	0	0	1	0	0	0	0;
         0	0	0	0	0	0	0	0	a1	a2	0	0;
         0	0	0	0	0	0	0	0	a2	a1	0	0;
         0	0	0	0	0	0	0	0	0	0	a1	a2;
         0	0	0	0	0	0	0	0	0	0	a2	a1];          

A =     [1	1	1	1	0	0	0	0	0	0	0	0;
         0	0	0	0	1	0	0	0	0	0	0	0;
         0	0	0	0	0	1	0	0	0	0	0	0;
         0	0	0	0	0	0	1	0	0	0	0	0;
         0	0	0	0	0	0	0	1	0	0	0	0;
         0	0	0	0	0	0	0	0	a3	a4	0	0;
         0	0	0	0	0	0	0	0	a4	a3	0	0;
         0	0	0	0	0	0	0	0	0	0	a3	a4;
         0	0	0	0	0	0	0	0	0	0	a4	a3];  

% set the strength of the uniform prior for A
a = ups*A;  

% The higher the multiplier, the slower the learning as priors are adjusted
% by numbers <1 on each trial. A very high order number would correspond to
% a very well learnt prior.

% B - (controlled) transitions between hidden states
%--------------------------------------------------------------------------
% Most transitions are deterministic, apart from 20% chance of aversive
% state even with the correct go/no-go response if learning in B is
% assumed. The last 4 states are absorbing. 
%--------------------------------------------------------------------------

% If learning in B is assumed or not
if strcmp(learnwhat,'b')
    b1 = 0.8; b2 = 1 - b1;   % process likelihoods 
    b3 = 0.5; b4 = 1 - b3;   % initial model likelihoods
else
    b1 = 1; b2 = 1 - b1;     % process likelihoods 
    b3 = 1; b4 = 1 - b3;     % initial model likelihoods
end


B_ENV{1} = [1	0	0	0	0	0	0	0	0	0	0	0;
            0	1	0	0	0	0	0	0	0	0	0	0;
            0	0	1	0	0	0	0	0	0	0	0	0;
            0	0	0	1	0	0	0	0	0	0	0	0;
            0	0	0	0	0	0   0	0	0	0	0	0;
            0	0	0	0	0	0	0	0	0	0	0	0;
            0	0	0	0	0	0	0	0	0	0	0	0;
            0	0	0	0	0	0	0	0	0	0	0	0;
            0	0	0	0	0	0	0	0	1	0	0	0;
            0	0	0	0	0	0	0	0	0	1	0	0;
            0	0	0	0	b2	b1	b2	b1	0	0	1	0;
            0	0	0	0	b1	b2	b1	b2	0	0	0	1];

% Go to location 2
B_ENV{2} = [0	0	0	0	0	0	0	0	0	0	0	0;
            0	0	0	0	0	0	0	0	0	0	0	0;
            0	0	0	0	0	0	0	0	0	0	0	0;
            0	0	0	0	0	0	0	0	0	0	0	0;
            1	0	0	0	1	0	0	0	0	0	0	0;
            0	1	0	0	0	1	0	0	0	0	0	0;
            0	0	1	0	0	0	1	0	0	0	0	0;
            0	0	0	1	0	0	0	1	0	0	0	0;
            0	0	0	0	0	0	0	0	1	0	0	0;
            0	0	0	0	0	0	0	0	0	1	0	0;
            0	0	0	0	0   0   0   0  	0	0	1	0;
            0	0	0	0	0   0   0	0	0	0	0	1];  

% Go to location 3
B_ENV{3} = [1	0	0	0	0	0	0	0	0	0	0	0;
            0	1	0	0	0	0	0	0	0	0	0	0;
            0	0	1	0	0	0	0	0	0	0	0	0;
            0	0	0	1	0	0	0	0	0	0	0	0;
            0	0	0	0	0	0	0	0	0	0	0	0;
            0	0	0	0	0	0	0	0	0	0	0	0;
            0	0	0	0	0	0	0	0	0	0	0	0;
            0	0	0	0	0	0	0	0	0	0	0	0;
            0	0	0	0	b1	b2	b1	b2	1	0	0	0;
            0	0	0	0	b2	b1	b2	b1	0	1	0	0;
            0	0	0	0	0	0	0	0	0	0	1	0;
            0	0	0	0	0	0	0	0	0	0	0	1];  

B = B_ENV;    

B{1}(11:12,5:8) = [b4 b3 b4 b3; 
                   b3 b4 b3 b4];

B{3}(9:10,5:8)  = [b3 b4 b3 b4;
                   b4 b3 b4 b3];        

 
for i = 1:numel(B)
    b{i} = ups.*B{i};
end

% B0 - set tranistion probabilities for habitual/Pavlovian responses
%--------------------------------------------------------------------------
% Assumes vigorous pavlovian responses to proximal threats (aversive sound)
% and passive responses to neutral context (avoid condition)
%--------------------------------------------------------------------------

b0 = z; 

B0 =   [0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	0	0	0	0	0	0;
        1	0	0	0	0	0	0	0	0	0	0	0;
        0	1	0	0	0	0	0	0	0	0	0	0;
        0	0	1	0	0	0	0	0	0	0	0	0;
        0	0	0	1	0	0	0	0	0	0	0	0;
        0	0	0	0	0	0	b0  b0	1	0	0	0;
        0	0	0	0	0	0 1-b0 1-b0	0	1	0	0;
        0	0	0	0   b0 b0	0	0	0	0	1	0;
        0	0	0	0 1-b0 1-b0 0	0	0	0	0	1];
    
% C - preferences of outcomes (utility)
%--------------------------------------------------------------------------
% Silence is neutral (0) while aversive outcomes are dispreferred (-c)
%--------------------------------------------------------------------------
C  = [0 0 0 -c -c 0 -c 0 -c]';

% D - prior beliefs about the context at t=1
%--------------------------------------------------------------------------
% Assumes that each of the 4 contexts are equally likely at the beginning
% of a trial
%--------------------------------------------------------------------------
d = [1 1 1 1 0 0 0 0 0 0 0 0]';

% Allowable instrumental policies 
% Pavlovian policy gets added later in spm_MDP_VB_LC_EA.m
%--------------------------------------------------------------------------  
V = [2 2;  % t = 1 - 'move' to location 2 is compulsory    
     3 1]; % t = 2 - a choice between go (3) and no-go (1)

% MDP Structure - this will be used to generate arrays for multiple trials
%==========================================================================

mdp.A = A;                     % observation model
mdp.B = B;                     % transition probabilities
mdp.C = C;                     % preferred states
mdp.V = V;                     % allowable policies
mdp.e = e;                     % preferred policies
mdp.B0 = B0;                   % transitions for habitual/fast responses


% set concentration parameters if learning is assumed
if strcmp(learnwhat,'b')
    mdp.b = b;
end

if strcmp(learnwhat,'a')
    mdp.a = a;
end

mdp.d = d;                     % prior over contexts
mdp.s = 1;                     % initial state
mdp.A_ENV = A_ENV;             % real world (fixed) states-outcomes matrix 
mdp.B_ENV = B_ENV;             % real world (fixed) state transition matrix
mdp.alpha = mdp.alpha;         % precision of action selection
mdp.beta  = mdp.beta;          % inverse precision of policy selection
mdp.Ni = 16;                   % number of iterative updates

if df==0
    mdp.df_set=[];
else
    mdp.df_set=df;
end

% store all other params for plotting purposes
mdp.c   = c;
mdp.k   = k;
mdp.df  = df;
mdp.ups = ups;

mdp.dg  = mdp.dg; 
mdp.dh  = mdp.dh;
mdp.dl  = mdp.dl;
mdp.dm  = mdp.dm;
mdp.eta = mdp.eta;
mdp.z   = mdp.z;
mdp.w0  = mdp.w0;

% generate trials
n   = n_trials;                           % number of trials
ff  = randperm(n);                        % permutate trial sequence
Av  = ff(1:round(n*(p_av)));              % Avoid trials
Es  = ff(round(n*(p_av))+1:end);          % Escape trials
GA  = Av(1:round(numel(Av)*(p_Go)));      % Go-Avoid trials;     s=1
NGA = Av(round(numel(Av)*(p_Go))+1:end);  % No-Go-Avoid trials;  s=2
GE  = Es(1:round(numel(Es)*(p_Go)));      % Go-Escape trials;    s=3
NGE = Es(round(numel(Es)*(p_Go))+1:end);  % No-Go-Escape trials; s=4
                               
mdp.rev_n = rev_n;               % store index of trial when probs reverse

[MDP(1:n)] = deal(mdp);          % sets MDP up with an expanded struct

% assign initial states to different contexts
for i = 1:n    
    if ismember(i, GA)        
        [MDP(i).s] = 1;
    elseif ismember(i, NGA) 
        [MDP(i).s] = 2;
    elseif ismember(i, GE) 
        [MDP(i).s] = 3;
    elseif ismember(i, NGE) 
        [MDP(i).s] = 4;
    end
end

% if a reversal happens, change the environment 
for i = rev_n+1:n    
       
    if strcmp(learnwhat,'b')
        MDP(i).B_ENV{1}(11:12,5:8) = [b1 b2 b1 b2; 
                                      b2 b1 b2 b1]; 

        MDP(i).B_ENV{3}(9:10,5:8) = [b2 b1 b2 b1; 
                                     b1 b2 b1 b2]; 
                                     
    elseif strcmp(learnwhat,'a')

        ap = a1; a1 = a2; a2 = ap; 

        MDP(i).A_ENV = [1	1	1	1	0	0	0	0	0	0	0	0;
                        0	0	0	0	1	0	0	0	0	0	0	0;
                        0	0	0	0	0	1	0	0	0	0	0	0;
                        0	0	0	0	0	0	1	0	0	0	0	0;
                        0	0	0	0	0	0	0	1	0	0	0	0;
                        0	0	0	0	0	0	0	0	a1	a2	0	0;
                        0	0	0	0	0	0	0	0	a2	a1	0	0;
                        0	0	0	0	0	0	0	0	0	0	a1	a2;
                        0	0	0	0	0	0	0	0	0	0	a2	a1];          

    end
        
end
    
end