function plot_GoNoGo_EA_many_many(mdp,opts)
% plots performance while varying a set of parameters (mdp.pv)
% one parameter at a time

os = mdp;           % save the original settings 
np = numel(mdp.pn); % number of parameters to vary
reps = mdp.runs;    % number of repetitions with the same params

%loop over parameters
for j = 1:np
    
    mdp = os;       % reset the param values after previous manipulation

    % extract summary statistics and plot full performance for each s.pv
    for i = 1:size(mdp.pv,2)

        % set the desired parameter
        mdp.(mdp.pn{j}) = mdp.pv(j,i);

        for ii = 1:reps

            MDP = GoNoGo_EA(mdp,mdp.learnwhat); % generate task structure
            rng('shuffle')                      % pick a random seed 
            MDP = spm_MDP_VB_LC_EA(MDP);        % simulate task performance 
            
            ss = plot_GoNoGo_EA_stats(MDP,opts);    % plot and get stats

            Acc_br(j,ii,i,:)   = ss.Acc_br;
            Acc_ar(j,ii,i,:)   = ss.Acc_ar;
            Acc(j,ii,i,:)      = ss.Acc;
            AEB(j,ii,i,:)      = ss.AEB;     
            ABNM(j,ii,i,:)     = ss.ABNM;
            ABGM(j,ii,i,:)     = ss.ABGM;
            RPM(j,ii,i,:)      = ss.RPM;
        end
    end
    
    % mean accuracy across conditions
    MAcc(j,:,:,:) = cat(3,squeeze(mean(Acc_br(j,:,:,:),4)),...
        squeeze(mean(Acc_ar(j,:,:,:),4)),squeeze(mean(Acc(j,:,:,:),4))); 
        
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
end

% Plot performance as a function of s.pv (before/after reversal and overall)
%--------------------------------------------------------------------------

figure('WindowStyle','docked');
%loop over parameters
for j = 1:np    
    
    rs = 4; % figure rows
    
    cols2 = [0.08,0.70,0.55;
             0.65,0.65,0.05]; % colors for before/after reversal
         
    subplot(rs,np,j); hold on
    plotbat_rep(mdp.pv(j,:),squeeze((ABGM(j,:,:,:)+ABNM(j,:,:,:))/2),cols2)
    ylim([0.2 0.6]);
    if j==1, ylabel('Mean beliefs'), end
    set(gca, 'FontSize',11)
    
    subplot(rs,np,j+np); hold on
    plotbat_rep(mdp.pv(j,:),squeeze(RPM(j,:,:,:)),cols2)
    ylim([0 0.37]); 
    if j==1, ylabel('Pavlovian policy'), end  
    set(gca, 'FontSize',11)

    subplot(rs,np,j+2*np); hold on
    ylim([-0.01 0.51]); 
    xlim([mdp.pv(j,1)   - 0.065*range(mdp.pv(j,:)) ...
          mdp.pv(j,end) + 0.065*range(mdp.pv(j,:))])
    h1 = hline(0.11,'k-'); h2 = hline(0.22,'k:');
    set(h1,'LineWidth',2,'Color',[0.77,0.02,0.16]);
    set(h2,'LineWidth',2,'Color',[0.77,0.02,0.16]);
    plotbat_rep(mdp.pv(j,:),squeeze(AEB(j,:,:,:)),cols2)
    if j==1, ylabel('Escape bias'), end
    set(gca, 'FontSize',11)
    
    subplot(rs,np,j+3*np); hold on
    plotbat_rep(mdp.pv(j,:),squeeze(MAcc(j,:,:,:)),cols2)
    ylim([0.4 0.9]); hline(0.5,'k--');
    if j==1, ylabel('Mean accuracy'), end
    set(gca, 'FontSize',11)
    xlabel(mdp.pn(j));
    
%     plotbat(s.pv(j,:),PAB,cols2)
%     ylim([-0.2 1]);
%     if j==1, ylabel('Avoid bias'), end
%     xlabel(s.pn(j));
%     hline(0,'k--');
end

subplot(rs,np,1)
legend('Before reversal', 'After reversal', 'Overall')


% same as above but plot variance instead of means
% figure('WindowStyle','docked');
% %loop over parameters
% for j = 1:np    
%     
%     rs = 4; % figure rows
%     
%     cols2 = [0.08,0.70,0.55;
%              0.65,0.65,0.05]; % colors for before/after reversal
%          
%         subplot(rs,np,j); hold on
%     plotbat_var_rep(s.pv(j,:),squeeze((ABGM(j,:,:,:)+ABNM(j,:,:,:))/2),cols2)
%     %ylim([0.2 0.6]);
%     if j==1, ylabel('Mean beliefs'), end
%     set(gca, 'FontSize',11)
%     
%     subplot(rs,np,j+np); hold on
%     plotbat_var_rep(s.pv(j,:),squeeze(RPM(j,:,:,:)),cols2)
%     %ylim([0 0.37]); 
%     if j==1, ylabel('Pavlovian policy'), end  
%     set(gca, 'FontSize',11)
% 
%     subplot(rs,np,j+2*np); hold on
%     %ylim([-0.01 0.51]); 
%     plotbat_var_rep(s.pv(j,:),squeeze(AEB(j,:,:,:)),cols2)
%     if j==1, ylabel('Escape bias'), end
%     set(gca, 'FontSize',11)
%     
%     subplot(rs,np,j+3*np); hold on
%     plotbat_var_rep(s.pv(j,:),squeeze(MAcc(j,:,:,:)),cols2)
%     %ylim([0.4 0.9]); hline(0.5,'k--');
%     if j==1, ylabel('Mean accuracy'), end
%     set(gca, 'FontSize',11)
%     xlabel(s.pn(j));
%     
% %     plotbat(s.pv(j,:),PAB,cols2)
% %     ylim([-0.2 1]);
% %     if j==1, ylabel('Avoid bias'), end
% %     xlabel(s.pn(j));
% %     hline(0,'k--');
% end
% 
% subplot(rs,np,1)
% legend('Before reversal', 'After reversal', 'Overall')



% tight subplots
% subtightplot(rs,np,j); hold on
% plotbat(s.pv(j,:),ABGM,cols2)
% ylim([0.1 0.7]);
% if j==1
%     ylabel('Mean of beliefs'), 
%     set(gca, 'XTickLabel', []);
% else
%     set(gca, 'XTickLabel', [],'YTickLabel', []);
% end
% 
% subtightplot(rs,np,j+np); hold on
% plotbat(s.pv(j,:),RPM,cols2)
% ylim([0 1]); 
% if j==1
%     ylabel('Mean Pavlovian policy')
%     set(gca, 'XTickLabel', []);
% else
%     set(gca, 'XTickLabel', [],'YTickLabel', []);
% end    
% 
% subtightplot(rs,np,j+2*np); hold on
% plotbat(s.pv(j,:),AEB,cols2)
% ylim([-0.2 1]); hline(0,'k--'); 
% if j==1
%     ylabel('Escape bias')
%     set(gca, 'XTickLabel', []);
% else
%     set(gca, 'XTickLabel', [],'YTickLabel', []);
% end    
% 
% subtightplot(rs,np,j+3*np); hold on
% plotbat(s.pv(j,:),PAB,cols2)
% ylim([-0.2 1]); hline(0,'k--');
% xlabel(s.pn(j));
% if j==1
%     ylabel('Avoid bias')
% else
%     set(gca,'YTickLabel', []);
% end    


% original params
c   = os.c;
k   = os.k;
al  = os.alpha;
bt  = os.beta;
e   = os.e;
ups = os.ups;
z   = os.z;
w0  = os.w0;

if isempty(MDP(1).df_set)
    df = nan;
else
    df = os.df;
end

dm = os.dm;
dg = os.dg;
dl = os.dl;
dh = os.dh;


% list all of the params in the suptitle
try 
    str_val=[];
    for i = 1:numel(opts.showparams)
        str_val = [str_val opts.showparams{i}, ' = ', string(eval(opts.showparams{i})), '   '];
    end
    sp = suptitle(sprintf('%s',str_val{1:end}));
catch
    sp = suptitle(sprintf('c = %g, k = %g, \\eta = %g, \\eta_g = %.2f, \\eta_m = %.2f, \\eta_{max} = %.2f, \\eta_{min} = %.2f, \\alpha = %g, \\beta = %.1f, e = [%g %g], ups = %g',...
        c, k, df, dg, dm, dh, dl, al, bt, e, ups));
end
set(sp,'FontWeight','Bold');

end

function plotbat_rep(X,Y,c)
% plots before and after reversal performance + overall performance

if ndims(Y)>2 %#ok<ISMAT>           % if more than one repeated run
    YY = squeeze(mean(Y,1));
    err = squeeze(std(Y,1)/sqrt(size(Y,1)));
    errorbar(X,YY(:,1),err(:,1),'*-','LineWidth',2,'color',c(1,:))
    errorbar(X,YY(:,2),err(:,2),'*-','LineWidth',2,'color',c(2,:))
    errorbar(X,YY(:,3),err(:,3),'k*-','LineWidth',2)
else                                % if 1 run    
    plot(X,Y(:,1),'*-','LineWidth',2,'color',c(1,:))
    plot(X,Y(:,2),'*-','LineWidth',2,'color',c(2,:))
    plot(X,Y(:,3),'k*-','LineWidth',2)
end        
xlim([X(1)-X(end)*0.05 X(end)*1.05]);
end

function plotbat_var_rep(X,Y,c)
% plots before and after reversal performance + overall performance
YY = squeeze((std(Y,1)).^2);

plot(X,YY(:,1),'*-','LineWidth',2,'color',c(1,:))
plot(X,YY(:,2),'*-','LineWidth',2,'color',c(2,:))
plot(X,YY(:,3),'k*-','LineWidth',2)
xlim([X(1)-X(end)*0.05 X(end)*1.05]);
end