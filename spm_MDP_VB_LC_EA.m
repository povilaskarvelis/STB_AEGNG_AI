function [MDP] = spm_MDP_VB_LC_EA(MDP)
% active inference and learning using variational Bayes (modified)

% A modified version of spm_MDP_VB (see https://www.fil.ion.ucl.ac.uk/spm/)
% and spm_MDP_VB_LC (see https://github.com/AnnaCSales/ActiveInference).

% Original code spm_MDP_VB - Copyright (C) 2005 Wellcome Trust Centre for
% Neuroimaging - Karl Friston; 
% Amended code spm_MDP_VB_LC - Anna Sales, 2018, University of Bristol.
% Current code spm_MDP_VB_LC_EA - Povilas Karvelis, 2022, CAMH, Toronto

% FORMAT [MDP] = spm_MDP_VB(MDP)
%
% MDP.S(N,1)        - true initial state
% MDP.V(T - 1,P)    - P allowable policies (control sequences)
%
% MDP.A(O,N)        - likelihood of O outcomes given N hidden states
% MDP.B{M}(N,N)     - transition probabilities among hidden states (priors)
% MDP.C(N,1)        - prior preferences   (prior over future outcomes)
% MDP.D(N,1)        - prior probabilities (prior over initial states)
%
% MDP.a(O,N)        - concentration parameters for A
% MDP.b{M}(N,N)     - concentration parameters for B
% MDP.c(N,N)        - concentration parameters for habitual B
% MDP.d(N,1)        - concentration parameters for D
% MDP.e(P,1)        - concentration parameters for u
%
% optional:
% MDP.s(1,T)        - vector of true states
% MDP.o(1,T)        - vector of observations 
% MDP.u(1,T)        - vector of actions
% MDP.w(1,T)        - vector of precisions
%
% MDP.alpha         - upper bound on precision (Gamma hyperprior - shape [1])
% MDP.beta          - precision over precision (Gamma hyperprior - rate  [1])
%
% produces:
%
% MDP.P(M,T)        - probability of emitting action 1,...,M at time 1,...,T
% MDP.Q(N,T)        - an array of conditional (posterior) expectations over
%                         N hidden states and time 1,...,T
% MDP.X             - and Bayesian model averages over policies
% MDP.R             - conditional expectations over policies
%
% MDP.un            - simulated neuronal encoding of hidden states
% MDP.xn            - simulated neuronal encoding of policies
% MDP.wn            - simulated neuronal encoding of precision (tonic)
% MDP.dn            - simulated dopamine responses (phasic)
% MDP.rt            - simulated reaction times
%==========================================================================
 
 
% if there are multiple trials ensure that parameters are updated
%--------------------------------------------------------------------------
if length(MDP) > 1
    
    for i = 1:length(MDP)
        
        % update concentration parameters
        %------------------------------------------------------------------
        if i > 1
            try  MDP(i).a = OUT(i - 1).a; end
            try  MDP(i).b = OUT(i - 1).b; end
            try  MDP(i).c = OUT(i - 1).c; end
            try  MDP(i).d = OUT(i - 1).d; end
            try  MDP(i).e = OUT(i - 1).e; end
            try  MDP(i).beta = OUT(i - 1).beta; end
        end
        
        % solve this trial
        %------------------------------------------------------------------
        OUT(i) = spm_MDP_VB_LC_EA(MDP(i));

    end
    MDP = OUT;    
    return
end
 
% set up and preliminaries
%==========================================================================
V   = MDP.V;                        % allowable policies (T - 1,Np)

% numbers of transitions, policies and states
%--------------------------------------------------------------------------
T   = size(MDP.V,1) + 1;            % number of transitions
Np  = size(MDP.V,2);                % number of allowable policies
Ns  = size(MDP.B{1},1);             % number of hidden states
Nu  = size(MDP.B,2);                % number of hidden controls
Nh  = size(MDP.e,1);                % index of habit
p0  = exp(-6);                      % smallest probability
q0  = 1/16;                         % smallest probability
 
% parameters of generative model and policies
%==========================================================================
 
% likelihood model (for a partially observed MDP implicit in G)
%--------------------------------------------------------------------------
A      = spm_norm(MDP.A + p0);   % model likelihood
A_ENV  = spm_norm(MDP.A_ENV);    % process likelihood
No     = size(A,1);              % number of outcomes
 
% parameters (concentration parameters): A
%--------------------------------------------------------------------------
if isfield(MDP,'a')
    qA = MDP.a + q0;                          
    qA = psi(qA) - ones(No,1)*psi(sum(qA));
    qA = log(spm_softmax(qA));              
    A = spm_softmax(qA);                   
else
    
    qA = log(spm_norm(A));
end
 
% transition probabilities (priors)
%--------------------------------------------------------------------------
for i = 1:Nu
    
    B{i} = spm_norm(MDP.B{i} + p0);     % model state transitions
    B_ENV{i} = spm_norm(MDP.B_ENV{i});  % process state transitions
        
    % parameters (concentration parameters): B
    %----------------------------------------------------------------------
    if isfield(MDP,'b')
        b     = MDP.b{i} + q0;
        sB{i} = spm_norm(b );
        rB{i} = spm_norm(b');
        qB{i} = psi(b) - ones(Ns,1)*psi(sum(b)); 
        B{i} = spm_softmax(qB{i});        
    else
        b     = MDP.B{i} + p0;
        sB{i} = spm_norm(b );
        rB{i} = spm_norm(b');
    end
end 
 
% priors over initial hidden states - concentration parameters
%--------------------------------------------------------------------------
if isfield(MDP,'d')
    d  = MDP.d + q0;                        
    qD = psi(d) - ones(Ns,1)*psi(sum(d));  
    qD = log(spm_softmax(qD));             
elseif isfield(MDP,'D')
    d  = MDP.D + q0;
    qD = log(spm_norm(d));
else
    d  = ones(Ns,1);
    qD = psi(d) - ones(Ns,1)*psi(sum(d));
    qD = log(spm_softmax(qD));             
end

% priors over policies - concentration parameters
%--------------------------------------------------------------------------
e   = MDP.e;
qE  = psi(e) - ones(Nh,1)*psi(sum(e));     
qE  = log(spm_softmax(qE));                 

% habitual/Pavlovian state transition matrix
%--------------------------------------------------------------------------
B0    = MDP.B0 + p0;
sH    = spm_norm(B0 );
rH    = spm_norm(B0');
%qH    = psi(B0) - ones(Ns,1)*psi(sum(B0));  % ln(B0) = psi(B0) - psi(B0_0)

% prior preferences (log probabilities) : C
%--------------------------------------------------------------------------
Vo = MDP.C;    

% assume constant preferences, if only final states are specified
%--------------------------------------------------------------------------
if size(Vo,2) ~= T
    Vo = Vo(:,end)*ones(1,T);
end

Vo    = log(spm_softmax(Vo)+eps); % log preferences ln(P(o))
H     = sum(spm_softmax(qA).*qA); % entropy 

% log preferences over states
%--------------------------------------------------------------------------
Vs    = spm_norm(A')*spm_softmax(Vo);
Vs    = log(Vs) + H'*ones(1,T);
Vs    = log(spm_softmax(Vs));
    
% precision defaults
%--------------------------------------------------------------------------
alpha = MDP.alpha;  
beta  = MDP.beta;
eta   = MDP.eta;
 
% initial states and outcomes
%--------------------------------------------------------------------------
try
    s = MDP.s(1);                   % initial state (index)
catch
    s = 1;
end

try
    o = MDP.o(1);                   % initial outcome (index)
catch
    o = find(rand < cumsum(A_ENV(:,s)),1);
end

P  = zeros(Nu,T - 1);               % posterior beliefs about control
x  = zeros(Ns,T,Nh) + 1/Ns;         % expectations of hidden states | policy
X  = zeros(Ns,T);                   % expectations of hidden states
u  = zeros(Nh,T);                   % expectations of policy
a  = zeros(1, T - 1);               % action (index)
 
% initialise priors over states
%--------------------------------------------------------------------------
for k = 1:Nh
    x(:,1,k) = spm_softmax(qD);
end
 
% expected rate parameter
%--------------------------------------------------------------------------
qbeta = beta;                       % initialise rate parameters
gu    = zeros(1,T)  + 1/qbeta;      % posterior precision (policy)

% solve
%==========================================================================
Ni    = MDP.Ni;                     % number of VB iterations
xn    = zeros(Ni,Ns,T,T,Np) + 1/Ns; % history of state updates
p     = 1:Nh;                       % allowable policies

for t = 1:T
        
    % reset states
    %----------------------------------------------------------------------        
    x = spm_softmax(log(x)/4); % why 4?
    
    % Variational updates (hidden states) under sequential policies
    %======================================================================
    F = zeros(Nh,T);
    for k = p
        
        for i = 1:Ni
            px    = x(:,:,k);            
            
            for j = 1:T
                
                % current state
                %----------------------------------------------------------
                qx   = log(x(:,j,k));                
                   
                % transition probabilities (without attention)
                %------------------------------------------------------
                if k > Np
                    fB  = sH;
                    bB  = rH;
                else
                    if j > 1, fB = sB{V(j - 1,k)}; end
                    if j < T, bB = rB{V(j    ,k)}; end
                end
                
                % evaluate free energy and gradients (v = dFdx)
                %----------------------------------------------------------
                v    = qx;
                if j <= t, v = v - qA(o(j),:)';           end
                if j == 1, v = v - qD;                    end
                if j >  1, v = v - log(fB*x(:,j - 1,k));  end

                % free energy and belief updating
                %----------------------------------------------------------
                F(k,j)  = -x(:,j,k)'*v;    
                
                if j <  T, v = v - log(bB*x(:,j + 1,k));  end
                
                px(:,j) = spm_softmax(qx - v/8);   % why 8?               
            end
            
            % hidden state updates
            %--------------------------------------------------------------
            x(:,:,k) = px;            
        end
    end
    
    % (negative path integral of) free energy of policies (Q)
    %======================================================================
    Q     = zeros(Nh,T);
    for k = p
        for j = t+1:T
            qx     = A*x(:,j,k);
            Q(k,j) = qx'*(Vo(:,j) - log(qx)) + H*x(:,j,k);
        end
    end    
    
    % variational updates - policies and precision
    %======================================================================
    F     = sum(F(:,1:t),2);
    Q     = sum(Q,2);

    % eliminate unlikely policies
    %----------------------------------------------------------------------
    p     = p(F(p) - max(F(p)) > -3);
    
    for i = 1:Ni        
        % policy (u)
        %------------------------------------------------------------------
        qu = spm_softmax(qE(p) + gu(t)*Q(p) + F(p)); 
        pu = spm_softmax(qE(p) + gu(t)*Q(p));        
        
        v     = qbeta - beta + (qu - pu)'*Q(p); 
        qbeta = qbeta - v/2;
        gu(t) = 1/qbeta;       
    end
    u(p,t)  = qu; % store policy probabilities (pi)
    
    % Bayesian model averaging of hidden states over policies
    %----------------------------------------------------------------------
    for i = 1:T
        X(:,i)     = squeeze(x(:,i,:))*u(:,t);
        X_t(:,i,t) = X(:,i);        
    end
    
    % SAPE = KL divergence between successive BMA distributions
    %----------------------------------------------------------------------
    if t > 1
        St_diff   = log(X_t(:,3,t)) - log(X_t(:,3,t-1)); 
        SAPE(t-1) = sum(sum(X_t(:,3,t).*St_diff));   
    end   
    
    % action selection and sampling of next state (outcome)
    %======================================================================
    if t < T
    
        % posterior expectations about (remaining) actions (q)
        %==================================================================
        q = 1:Nu;
        
        v     = log(A*X(:,t + 1));          
        for j = q
            qo     = A*B{j}*X(:,t);          
            P(j,t) = (v - log(qo))'*qo + 16; % why 16? 
        end
                
        % action selection
        %------------------------------------------------------------------        
        
        % a correction due to added smallest probabilities
        if t==1
            P(2,t) = exp(10);  % at t=1 only B(2) is possible
        elseif t==2
            P(2,t) = 0;        % at t=2 B(2) is not possible
        end

        % probability matching for action selection
        P(:,t) = spm_softmax(alpha*P(:,t));       
        a(t)   = find(rand < cumsum(P(:,t)),1);
            
        % next sampled state
        %------------------------------------------------------------------
        s(t + 1) = find(rand < cumsum(B_ENV{a(t)}(:,s(t))),1);
       
        % next observed state
        %------------------------------------------------------------------
        o(t + 1) = find(rand < cumsum(A_ENV(:,s(t + 1))),1);
       
        % gamma
        %------------------------------------------------------------------
        gu(1,t + 1)   = gu(t);
                
    end    
end

% learning
%==========================================================================
if isfield(MDP, 'df_set') && length(MDP.df_set)==1 % fixed df  
    df = MDP.df_set;
    SAPE_df = SAPE(2);
else                                               % SAPE-dependent df 
    dm      = MDP.dm; 
    min_d   = MDP.dl; 
    max_d   = MDP.dh; 
    grad_d  = MDP.dg;
    SAPE_df = SAPE(2);
    df      = logist(SAPE_df,grad_d, max_d, min_d, dm); 
end

% controllability 
eo     = A*X_t(:,3,2);                           % expected outcomes at tau=2, for t=3 
w      = logist(sum(eo([6,8])),-10,1,0,MDP.w0);  % controllability
MDP.CP = MDP.C.*(1-w);                           % modulate stress sensitivity

for t = 1:T
    
    % mapping from hidden states to outcomes: a
    %----------------------------------------------------------------------
    if isfield(MDP,'a')               
        i          = MDP.a > 0; 
        da         = sparse(o(t),1,1,No,1)*X(:,t)';
        ii         = false(size(MDP.a)); 
        ii(o(t),:) = i(o(t),:); 
        MDP.a(ii)  = MDP.a(ii) - ((MDP.a(ii) - 1)/df);   % decay
        MDP.a(ii)  = MDP.a(ii) + da(ii)*eta;             % update 
    end
    
    % mapping from hidden states to hidden states: b(u)
    %----------------------------------------------------------------------
    if isfield(MDP,'b') && t > 2   
        
        % learning rate with a potential boost
        et = eta + MDP.k*abs(MDP.CP(o(t))); 
        
        % Optional: policy blending (Pavlovian + instrumental)
%         uc = u;        
%         if any(s(2)==[5,6])
%            uc(1,2) = u(1,2) + u(3,2);
%         elseif any(s(2)==[7,8])
%            uc(2,2) = u(2,2) + u(3,2);
%         end
%         uc(3,2) = 0;        
%         u = uc;
            
        % updates 
        for k = 1:Np
            v            = V(t - 1,k);
            db           = u(k,t - 1)*x(:,t,k)*x(:,t - 1,k)';
            ii           = false(size(MDP.b{v})); 
            ii([9,10]+(k-1)*2,s(2)) = [1 1]; 
            MDP.b{v}(ii)  = MDP.b{v}(ii) - (MDP.b{v}(ii)-1)/df;
            MDP.b{v}(ii)  = MDP.b{v}(ii) + et*db(ii);
        end
    end    
end

% initial hidden states:
%--------------------------------------------------------------------------
if isfield(MDP,'d')
    i        = MDP.d > 0;
    MDP.d(i) = MDP.d(i) + X(i,1) - (MDP.d(i) - 1)/16  ;  
end

% assemble results and place in MDP structure
%--------------------------------------------------------------------------
MDP.P   = P;              % probability of action at time 1,...,T - 1
MDP.Q   = x;              % conditional expectations over N hidden states
MDP.X   = X;              % Bayesian model averages
MDP.R   = u;              % conditional expectations over policies
MDP.o   = o;              % outcomes at 1,...,T
MDP.s   = s;              % states at 1,...,T
MDP.u   = a;              % action at 1,...,T
MDP.w   = gu;             % posterior expectations of precision (policy)
MDP.C   = Vo;             % utility
MDP.Vs  = Vs;             % log preferences over states
MDP.T   = T;              % trial depth 
MDP.Ni  = Ni;             % number of iterations for VB
 
MDP.A       = A;          % mapping from states to outcomes
MDP.B       = B;          % state transition probabilities
MDP.X_t     = X_t;        % BMA at each time point
MDP.SAPE_df = SAPE_df;    % SAPE used as input for decay factor    
MDP.df      = df;         % decay factor used in trial
MDP.beta    = 1/gu(T);    % carry forward beta
%MDP.rt      = rt;         % processing reaction time

% MDP.un  = un;             % simulated neuronal encoding of policies
% MDP.xn  = Xn;             % simulated neuronal encoding of hidden states
% MDP.wn  = wn;             % simulated neuronal encoding of precision
% MDP.dn  = dn;             % simulated dopamine responses (deconvolved)
end

function A = spm_norm(A)
% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
A  = A*diag(1./sum(A,1));
end

function [y] = logist(x, k, maxv, minv, meanv)  
% calculate the decay factor df
    maxv = maxv - minv;
    y = minv + (maxv./(1+exp(k*(x-meanv))));
end