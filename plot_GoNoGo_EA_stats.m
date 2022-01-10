function ss = plot_GoNoGo_EA_stats(MDP,opts)
% -------------------------------------------------------------------------
% Plots statistical breakdown of behavior in the Go/No-Go Escape/Avoid task

% Extract relevant variables
%--------------------------------------------------------------------------
for i = 1:length(MDP) 
    contexts(i,:) = MDP(i).s(1);          % 1-GA, 2-NGA, 3-GE, 4-NGE
    actions(i,:)  = MDP(i).u(end);        % 1 - No-go, 3 - Go   
    try
        beliefs_states(i,:) = exp(MDP(i).d(1:4));   % [GA NGA GE NGE]
    catch
        beliefs_states(i,:) = zeros(1,4);
    end
    action_probs(i,:) = MDP(i).P(:,2);     % [null No-go Go]
    outcomes(i,:)     = MDP(i).o(end);  
    policies(i,:)     = MDP(i).R(:,2);     % policy probabilities at t=2
    
    dft(i,:)      = MDP(i).df;             % decay factor
    SAPE_df(i,:)  = MDP(i).SAPE_df';       % SAPE (used for computing df)
    
    B_nogo(i,:) = [MDP(i).B{1}(11,5), MDP(i).B{1}(11,6), MDP(i).B{1}(11,7), MDP(i).B{1}(11,8)];
    B_go(i,:)   = [MDP(i).B{3}(9,5), MDP(i).B{3}(9,6), MDP(i).B{3}(9,7), MDP(i).B{3}(9,8)];
end

beliefs_states = beliefs_states./sum(beliefs_states,2); % normalize beliefs (D)
ED  = -sum(beliefs_states.*log(beliefs_states),2);      % entropy of beliefs (D)

ABG = mean(B_go,2);
ABN = mean(B_nogo,2);

n_trials = length(MDP);       % the number of trials 
rev_n    = MDP(1).rev_n;      % reversal trial index
V        = MDP(1).V;          % allowable policies

% unpack parameters
c   = MDP(1).c;
k   = MDP(1).k;
al  = MDP(1).alpha;
bt  = MDP(1).beta;
e   = MDP(1).e;
ups = MDP(1).ups;
z   = MDP(1).z;
w0  = MDP(1).w0;

if isempty(MDP(1).df_set)
    df = nan;
else
    df = MDP(1).df;
end

dm = MDP(1).dm;
dg = MDP(1).dg;
dl = MDP(1).dl;
dh = MDP(1).dh;

% Compute various summary statistics
%--------------------------------------------------------------------------
% Contexts: GA - 1, NGA - 2, GE - 3, NGE - 4
% Actions: Go - 3, No-go - 2;

cab = [V(2,:),V(2,:)]; % correct actions before reversal (GA, NGA, GE, NGE)
caa = flip(cab);       % correct actions after reversal (GA, NGA, GE, NGE)
    
% in case there is no reversal
if rev_n > n_trials  
    for j = 1:4  % loop over GA, NGA, GE, NGE
        Acc_br(j) = mean(actions(contexts==j)==cab(j)); % accuracy
    end  
    EDM = mean(ED);                        % mean entropy 
    EDV = var(ED);                         % variance entropy
else
    cb = contexts(1:rev_n);                % contexts before reversal
    ca = contexts(rev_n+1:end);            % contexts after reversal
    ab = actions(1:rev_n);                 % actions before reversal
    aa = actions(rev_n+1:end);             % actions after reversal
    
    % loop over GA, NGA, GE, NGE
    for j = 1:4 
        Acc_br(j) = mean(ab(cb==j)==cab(j)); % accuracy before reversal
        Acc_ar(j) = mean(aa(ca==j)==caa(j)); % accuracy after reversal
        
        % Overal accuracy and rt in the task
        Acc(j) = (Acc_br(j)*rev_n + Acc_ar(j)*(n_trials-rev_n))/n_trials;
    end
    
    % Get summary statistics across conditions
    [EDM, EDV] = sum_stats(ED,rev_n,n_trials);                % Belief entropy
    [SAPE_dfm, SAPE_dfv] = sum_stats(SAPE_df,rev_n,n_trials); % SAPE
    [dftm, dftv] = sum_stats(dft,rev_n,n_trials);             % Decay factor
    [ABGM, ABGV] = sum_stats(ABG,rev_n,n_trials);             % Transition beliefs
    [ABNM, ABNV] = sum_stats(ABN,rev_n,n_trials);             % Transition beliefs
    [RPM, RPV]   = sum_stats(policies(:,3),rev_n,n_trials);   % Reflexive policy 
    
    % flip actions after reversal for plotting 'correct action probability'
    action_probs(rev_n+1:end,flip(V(2,:))) = action_probs(rev_n+1:end,V(2,:));
end

% active escape bias
AEB(1,:) = Acc_br(:,3) - Acc_br(:,4);                             % before
AEB(2,:) = Acc_ar(:,4) - Acc_ar(:,3);                             % after
AEB(3,:) = (AEB(1,:)*rev_n + AEB(2,:)*(n_trials-rev_n))/n_trials; % overall

% passive avoid bias
PAB(1,:) = Acc_br(:,2) - Acc_br(:,1);                             % before
PAB(2,:) = Acc_ar(:,1) - Acc_ar(:,2);                             % after
PAB(3,:) = (PAB(1,:)*rev_n + PAB(2,:)*(n_trials-rev_n))/n_trials; % overall


% store summary statistics in the ouput structure
ss.SAPE_dfm = SAPE_dfm;
ss.SAPE_dfv = SAPE_dfv;
ss.Acc_br   = Acc_br;
ss.Acc_ar   = Acc_ar;
ss.Acc      = Acc;
ss.AEB      = AEB;
ss.PAB      = PAB;
ss.EDM      = EDM;
ss.EDV      = EDV;
ss.dftm     = dftm;
ss.dftv     = dftv;
ss.ABNM     = ABNM;
ss.ABGM     = ABGM;
ss.ABNV     = ABNV;
ss.ABGV     = ABGV;
ss.RPM      = RPM;
ss.RPV      = RPV;

% Plotting 
%--------------------------------------------------------------------------
ctx_names = {'GA','NGA','GE','NGE'};
if opts.mainplot == 1

    figure('Name', 'Mainplot','WindowStyle','docked');
    cols_def = get(gca,'colororder');
    cols  = cols_def([1,2,5,3],:);  % colors for GA NGA GE NGE
    cols2 = cols_def([6,7],:);      % colors for before/after reversal
    
    % Accuracy before reversal
    subplot(3,4,1); hold on
    xlim([0 5]); %hline(mean(Acc_br),'k--');
    daBarPlot(Acc_br,'groups',...
        [1,2,3,4],'colors',cols,'xtlabels',ctx_names);
    ylabel('Accuracy'); ylim([0 1]); yticks([0 0.2,0.4,0.6,0.8,1]);    

    if rev_n < n_trials
        title('Before reversal')
        % Accuracy after reversal
        subplot(3,4,2); hold on 
        xlim([0 5]); %hline(mean(Acc_ar),'k--');
        daBarPlot(Acc_ar,'groups',...
            [1,2,3,4],'colors',cols,'xtlabels',ctx_names(1,[2,1,4,3]));
        ylabel('Accuracy'); ylim([0 1]); yticks([0 0.2,0.4,0.6,0.8,1]);
        title('After reversal');
        
        % Overall performance accuracy
        subplot(3,4,3); hold on 
        xlim([0 5]); %hline(mean(Acc),'k--');
        daBarPlot(Acc,'groups',[1,2,3,4],...
            'colors',cols,'xtlabels',{' GA\newlineNGA',...
            'NGA\newline GA',' GE\newlineNGE','NGE\newline GE'});
        ylabel('Accuracy'); ylim([0 1]); yticks([0 0.2,0.4,0.6,0.8,1]);
        title('Overall');
    end    
    
    % SAPE vs learning rate
    avno_indx(1,:) = outcomes==6 | outcomes==8; % neutral outcome indices
    avno_indx(2,:) = outcomes==7 | outcomes==9; % aversive outcome indices
    
    subplot(3,4,4); hold on
    scatter(SAPE_df(avno_indx(1,:)),dft(avno_indx(1,:)),...
        'MarkerEdgeColor', cols2(1,:));
    scatter(SAPE_df(avno_indx(2,:)),dft(avno_indx(2,:)),...
        'MarkerEdgeColor', cols2(2,:));
    xlim([0.25 2.5]); ylim([0 dh*1.05]);
    ylabel('Decay parameter'); xlabel('SAPE');
    legend('Neutral outcome','Aversive outcome');
    title('Decay rate');    
    
    % Action probabilities over trials for Go and No-go in Avoid condition
    
    % prepare patch coordinates for showing outcomes of each trial
    X_patch = repmat([0.5; 0.5; 1.5; 1.5],1,n_trials) + repmat(0:n_trials-1,4,1);
    Y_patch = repmat([1.2; 1.3; 1.3; 1.2],1,n_trials);
    
    subplot(3,1,2); hold on
    for j = 1:size(beliefs_states,2)
        patch(X_patch(:,find(contexts==j)),Y_patch(:,find(contexts==j)),cols(j,:))
    end
    
    % add the final outcomes (aversive stimuli or no stimuli)        
    gav_indx(1,:) = outcomes==6; % go neutral outcomes
    gav_indx(2,:) = outcomes==7; % go aversive outcomes
    gav_indx(3,:) = outcomes==8; % nogo neutral outcomes
    gav_indx(4,:) = outcomes==9; % nogo aversive outcomes
    
    cols2 = [1, 1, 1;           % white (neutral outcome)
             0.1, 0.1, 0.1];    % black-dark grey (aversive outcome)
    
     patch(X_patch,Y_patch+0.1,[0.6 0.6 0.6]);
     patch(X_patch,Y_patch-0.1,[0.6 0.6 0.6]);
     h{1} = patch(X_patch(:,gav_indx(1,:)),Y_patch(:,gav_indx(1,:))+0.1,cols2(1,:));
     h{2} = patch(X_patch(:,gav_indx(2,:)),Y_patch(:,gav_indx(2,:))+0.1,cols2(2,:));
     patch(X_patch(:,gav_indx(3,:)),Y_patch(:,gav_indx(3,:))-0.1,cols2(1,:));
     patch(X_patch(:,gav_indx(4,:)),Y_patch(:,gav_indx(4,:))-0.1,cols2(2,:));
    
    for j = 1:4 % loop over GA, NGA, GE, NGE
        plot(find(contexts==j),action_probs(find(contexts==j),cab(j)),...      
            'color',cols(j,:), 'LineWidth',2)
    end     
    ylabel('Probability'); ylim([0 1]); yticks([0 0.2,0.4,0.6,0.8,1]);
    vline(rev_n,'k--'); xlim([0 n_trials+0.5]); ylim([0 1.4]);
    xlabel('Trials')
    
    % add scaled df value trajectory
    %h{3} = plot(1:n_trials,dft./max(dft),'k.-');    
    h{3} = plot(movmean(dft./max(dft),[3 3]),'k.-');
    
    title('Correct action probabilities'); 
        legend([h{:}],{'Neutral outcome','Aversive noutcome',...
        'Decay parameter \newline(scaled)'})

    % evolution of transition probs and policy probs
    ntr = 1:n_trials;
    for i = 1:4
     
        subplot(3,4,8+i); hold on
        ylim([0 1]);
        
        hb{1} = plot(ntr(contexts==i),policies(contexts==i,3),'color','k','LineWidth',2);
        hb{2} = plot(ntr(contexts==i),policies(contexts==i,2),'color',[0.5 0.5 0.5],'LineWidth',2,'Linestyle','-.');
        hb{3} = plot(ntr(contexts==i),policies(contexts==i,1),'color',[0.5 0.5 0.5],'LineWidth',2);
        hb{4+2*(i-1)} = plot(ntr(contexts==i),B_go(contexts==i,i),'color',cols(i,:),'LineWidth',2);
        hb{5+2*(i-1)} = plot(ntr(contexts==i),B_nogo(contexts==i,i),'color',cols(i,:),'LineWidth',2,'LineStyle','-.');

        xlim([0 n_trials]); xlabel('Trials'); vline(rev_n,'k--');
        cn = ctx_names(1,[2,1,4,3]); title([ctx_names{i} ' / ' cn{i}])
        
        if i==1               
            ylabel('Probability');
        elseif i==4
            legend([hb{4:end}, hb{3:-1:1}],...
                'B\{Go\}(s_9|s_5)','B\{NoGo\}(s_{11}|s_5)',...
                'B\{Go\}(s_9|s_6)','B\{NoGo\}(s_{11}|s_6)',...
                'B\{Go\}(s_9|s_7)','B\{NoGo\}(s_{11}|s_7)',...
                'B\{Go\}(s_9|s_8)','B\{NoGo\}(s_{11}|s_8)',...
                'Go policy','NoGo policy', 'Pavlovian policy',...
                'NumColumns',6);
        end       
    end 
    
    % list all of the params in the suptitle
    try 
        str_val=[];    
        for i = 1:numel(opts.showparams)
            str_val = [str_val, '   ',  opts.showparams{i}, ' = ', string(eval(opts.showparams{i}))];
        end
        
        sp = suptitle(sprintf('%s',str_val{1:end}));
    catch
        sp = suptitle(sprintf('c = %g, k = %g, \\eta = %g, \\eta_g = %.2f, \\eta_m = %.2f, \\eta_{max} = %.2f, \\eta_{min} = %.2f, \\alpha = %g, \\beta = %.1f, e = [%g %g], ups = %g',...
            c, k, df, dg, m, dh, dl, al, bt, e, ups));
    end
    set(sp,'FontWeight','Bold');
end

end

function [YM, YV] = sum_stats(X,rev_n,n_trials)
%--------------------------------------------------------------------------
% Get means and variances before and after reversal + overall

YM(1) = nanmean(X(1:rev_n));
YM(2) = nanmean(X(rev_n+1:end));
YV(1) = nanvar(X(1:rev_n));
YV(2) = nanvar(X(rev_n+1:end));

YM(3) = (YM(1)*rev_n + YM(2)*(n_trials-rev_n))/n_trials;
YV(3) = (YV(1)*rev_n + YV(2)*(n_trials-rev_n))/n_trials;

end