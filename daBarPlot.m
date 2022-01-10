function h = daBarPlot(Y,varargin)
%daBarPlot draws neat bars for multiple groups and multiple conditions 
%
% Description:
%
%   Creates bars organized by condition and colored by group. Does not
%   require the groups to be of the same size. Has some essential options
%   (see below) and exports the main handles for further customization. 
%
% Syntax:
%
%   daBarPlot(Y)
%   daBarPlot(Y,param,val,...)
%   h = daBarPlot(Y)
%   h = daBarPlot(Y,param,val,...)
%
% Input Arguments:
%
%   Y - data input (matrix or cell array) containing all conditions and all
%   groups. If Y is a matrix, each column has to correspond to different
%   condition, while the groups need to be specified in 'groups' vector.
%   If Y is a cell array, each cell has to contain data matrices for each 
%   group (columns being different conditions). In such case, the grouping 
%   is done automatically based on the cell structure.   
%
% Optional Input Parameter Name/Value Pairs:
%
%   NAME              VALUE
%
%   'groups'          A vector containing grouping variables. By default
%                     assumes a single group for a matrix data input. 
%
%   'fill'            0 - non-filled bars (contrours only)
%                     1 - bars filled with color (default)
%
%   'colors'          The RGB matrix for box colors of different groups
%                     (each row corresponding to a different group). If
%                     boxplots are filled, these are the fill colors with 
%                     the edges being black. If boxplots are not filled,
%                     these colors are used for edges. These colors can be 
%                     also used for scatter data instead (see 'flipcolors')
%                     Default colors: default matlab colors. 
%   
%   'errorbars'       Draws errorbars on each of the bars
%                     0     - no errorbars
%                     'SD'  - standard deviation (default) 
%                     'SEM' - standard error of the mean
%                     'CI'  - 95% confidence interval
%
%   'scatter'         0 - no datta scatter (deffault)
%                     1 - on top of the boxplot 
%
%   'scattersize'     Size of the scatter markers. Default: 15
%
%   'scattercolors'   Colors for the scattered data: {face, edge}
%                     Default: {'k','w'}
%
%   'flipcolors'      Will flip the colors of scatter points and bars
%                     0 - bars colored by group (default)
%                     1 - scatter is colored by group
%
%   'scatteralpha'    Transparency of scattered data (between 0 and 1)
%                     Default: 1 (completely non-transparent)
%
%   'jitter'          0 - do not jitter scattered data 
%                     1 - jitter scattered data (default)
% 
%   'outliers'        Highlights the outliers in the plot. The outliers 
%                     are values below Q1-1.5*IQR and above Q3+1.5*IQR.
%                     0 - do not highlight outliers  
%                     1 - highlight outliers (default)
%
%   'symbol'          Symbol and color for highlighting outliers.
%                     Default: 'rx' (red crosses).
%
%   'baralpha'        Bar transparency (between 0 and 1)
%                     Default: 1 (completely non-transparent)
%
%   'barspacing'      Scales spacing between bars in the same condition. 
%                     Note that spacing is also dependent on bar width
%                     Default: 1
%
%   'barwidth'        Scales the width of all bars
%                     Default: 1
%
%   'xtlabels'        Xtick labels (a cell of chars) for conditions. If
%                     there is only 1 condition and multiple groups, then 
%                     xticks and xtlabels will automatically mark different
%                     groups.
%                     Default: conditions/groups are numbered in the input 
%                     order
%
%   'legend'          Names of groups (a cell) for creating a legend
%                     Default: no legend
%
%
% Output Arguments:
%
%   h - a structure containing handles for further customization of the 
%   produced plot:
%       cpos - condition positions
%       gpos - group positions
%       
%       graphics objects:
%       br - bars
%       sc - scattered data markers
%       ot - outlier markers
%       er - error bars 
%       lg - legend
%
%
%
%
% Povilas Karvelis <karvelis.povilas@gmail.com>
% 15/04/2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

h = struct;
p = inputParser;

% specify default options
addOptional(p, 'groups', []);
addOptional(p, 'fill', 1); 
addOptional(p, 'colors', get(gca,'colororder'));
addOptional(p, 'errorbars', 'SD');
addOptional(p, 'scatter', 0); 
addOptional(p, 'scattersize', 15)
addOptional(p, 'scattercolors', {'k','w'}); 
addOptional(p, 'flipcolors', 0);
addOptional(p, 'scatteralpha', 1); 
addOptional(p, 'jitter', 1);
addOptional(p, 'mean', 1);
addOptional(p, 'outliers', 0); 
addOptional(p, 'symbol', 'rx'); 
addOptional(p, 'baralpha', 1);
addOptional(p, 'barspacing', 1);
addOptional(p, 'barwidth', 1);
addOptional(p, 'linkline',0);
addOptional(p, 'xtlabels', []);
addOptional(p, 'legend', []);


% parse the input options
parse(p, varargin{:});
confs = p.Results;    

% get group indices and labels
if ~isempty(confs.groups)
    [Gi,Gn,Gv] = grp2idx(confs.groups);
    num_groups = numel(Gv);
end

% find the number of groups
if iscell(Y)    
    num_groups = numel(Y);
    
    y = []; Gi = [];    
    for g = 1:num_groups
        y = [y; Y{g}];
        Gi = [Gi; g*ones(size(Y{g},1),1)];
    end   
       
    % default numbered group labels
    if ~exist('Gn','var')
        for g = 1:num_groups
            Gn{g} = num2str(g);
        end
    end    
    Y = y; % replace the cell with a data array
    
elseif ismatrix(Y)     
    % assume 1 group if none are specified
    if isempty(confs.groups)
       Gi = ones(size(Y,1),1);
       num_groups = 1;
    end
end

% find condition positions
if any(size(Y)==1)
    Y = Y(:);
    cpos = 1;
else
    cpos = 1:size(Y,2);
end
num_locs = numel(cpos);

% use condition positions to scale spacings
gpos=[];
if num_locs==1
    gpos = (1:num_groups)';
    bar_width = 4/5*confs.barwidth;
    cpos=gpos;
else    
    if num_groups==1
        gpos = cpos;
        bar_width = 4/5*confs.barwidth;
    else
        bar_width = (2/3)/(num_groups+1)*confs.barwidth;  % bar width 
        loc_sp = (bar_width/3)*confs.barspacing; % local spacing of bars

        % set group positions for each group
        for g = 1:num_groups
            gpos = [gpos; cpos + (g-(num_groups+1)/2)*(bar_width + loc_sp)];
        end
    end
end

h.gpos = gpos;
h.cpos = cpos; 

% loop over groups
for g = 1:num_groups 
    
    % get percentiles
    pt = prctile(Y(Gi==g,:),[2 9 25 50 75 91 98]); 
    means = mean(Y(Gi==g,:));
    medians = median(Y(Gi==g,:));
    SEMs = var(Y(Gi==g,:))/sqrt(size(Y(Gi==g,:),1));
    
    
    if size(pt,1)==1 pt=pt'; end % a fix for plotting one condition
    
    IQR = (pt(5,:)-pt(3,:));
        
    % create coordinates for drawing bars
    yb = reshape([zeros(1,length(means)); zeros(1,length(means))], 1, []);
    yt = reshape([means; means], 1, []);

    x1 = [gpos(g,:) - bar_width/2; gpos(g,:) - bar_width/2];
    x2 = [gpos(g,:) + bar_width/2; gpos(g,:) + bar_width/2];

    bar_ycor = [yt; yb];        
    bar_xcor = reshape([x1; x2],2,[]); 
    box_mdcor = reshape([pt(4,:); pt(4,:)], 1, []);
    box_mncor = reshape([means; means], 1, []);
    
    % create coordinates for drawing whiskers with cross-hatches and ends    
    hat_xcor = [gpos(g,:) - bar_width/4; gpos(g,:) + bar_width/4];    
    whi_xcor = [gpos(g,:); gpos(g,:)];        
    
    % draw one box at a time
    for k = 1:num_locs
        
        data_vals = Y(Gi==g,k); % data for a single bar
        
        % determine outliers and whisker length 
        ol = data_vals<(pt(3,k)-IQR(k)); % indices of lower outliers
        ou = data_vals>(pt(5,k)+IQR(k)); % indices of upper outliers    

%         whi_ycor(:,1,k) = [min(data_vals(~ol)), pt(3,k)]; % lower whisker        
%         whi_ycor(:,2,k) = [max(data_vals(~ou)), pt(5,k)]; % upper whisker
        
        whi_ycor(:,1,k) = [means(k), means(k)-SEMs(k)]; % lower whisker        
        whi_ycor(:,2,k) = [means(k), means(k)+SEMs(k)]; % upper whisker
        
        
        % jitter or not
        if confs.jitter==1
            xdata =  gpos(g,k).*ones(numel(Y(Gi==g,k)),1) + ...
                (bar_width/3).*(0.5 - rand(numel(Y(Gi==g,k)),1));
        elseif confs.jitter==0
            xdata = gpos(g,k).*ones(numel(Y(Gi==g,k)),1);
        end        

        % index values for each box
        wk = (1:2)+2*(k-1);
        Xx = bar_xcor(1:2,wk); 
        Yy = bar_ycor(1:2,wk); 

        % filled or not filled boxes
        if confs.fill==0            
            % no fill box
            h.bx(k,g) = line([Xx(:,1)' Xx(1,:) Xx(:,2)' Xx(2,:)],...
                [Yy(:,1)' Yy(1,:) Yy(:,2)' Yy(2,:)],...
                'color',confs.colors(g,:),'LineWidth',1.5); 
            hold on;        
            
%             % draw the median
%             h.md(k,g) = line(Xx(1,:), box_mdcor(wk),...
%                 'color',confs.colors(g,:), 'LineWidth', 2);            
%             
%             % draw the mean
%             if confs.mean==1
%                 h.mn(k,g) = line(Xx(1,:),box_mncor(wk),'LineStyle',':',...
%                     'color',confs.colors(g,:),'LineWidth', 1.5);
%             end           
            
        elseif confs.fill==1
            % box filled with color 
            h.bx(k,g) = fill([Xx(:,1)' Xx(1,:) Xx(:,2)' Xx(2,[2,1])],...
                 [Yy(:,1)' Yy(1,:) Yy(:,2)' Yy(2,:)],confs.colors(g,:));            
            set(h.bx(k,g),'FaceAlpha',confs.baralpha); 
            hold on;

%             % draw the median
%             h.md(k,g) = line(Xx(1,:), box_mdcor(wk),...
%                 'color','k', 'LineWidth', 2);
%             
%             % draw the mean
%             if confs.mean==1
%                 h.mn(k,g) = line(Xx(1,:),box_mncor(wk),'LineStyle',':',...
%                     'color','k','LineWidth', 1.5);
%             end
        end        
        
        ox = data_vals>max(data_vals); % default - no outliers
        
        % draw outliers
        if confs.outliers==1            
            ox = data_vals<whi_ycor(1,1,k) | data_vals>whi_ycor(1,2,k);
            h.ot(k,g) = scatter(xdata(ox),data_vals(ox),confs.scattersize,...
                confs.symbol);            
        end        
        
        % draw error bars
        if confs.errorbars==1
            
            h.wh(k,g,:) = plot(whi_xcor(:,k),whi_ycor(:,1,k),'k-',... 
                hat_xcor(:,k),[whi_ycor(1,1,k) whi_ycor(1,1,k)],'k-',... 
                whi_xcor(:,k),whi_ycor(:,2,k),'k-',... 
                hat_xcor(:,k),[whi_ycor(1,2,k) whi_ycor(1,2,k)],'k-',... 
                'LineWidth',1);                              
        end 

        % scatter on top of the boxplots
        if confs.scatter==1 || confs.scatter==2            
            h.sc(k,g) = scatter(xdata(~ox),data_vals(~ox),...
                confs.scattersize,...
                'MarkerFaceColor', confs.scattercolors{1},...
                'MarkerEdgeColor', confs.scattercolors{2},...
                'MarkerFaceAlpha', confs.scatteralpha); 
            hold on;    
        end        
        
    end        
        
    % put scattered data underneath boxplots
    if confs.scatter==2
        uistack(h.sc(:,g),'bottom')
    end   
    
    if confs.linkline==1
       h.ln(g) = line(gpos(g,:),pt(4,:),'color',confs.colors(g,:),...
           'LineStyle','-.','LineWidth',1.5); 
    end
    
end

% move lines to the background
if confs.linkline==1
    uistack(h.ln,'bottom')
end


% flip scatter and box colors and make a legend
if confs.flipcolors==1    
    
    box_class = class(h.bx); % box filled or no
    
    if strcmp(box_class,'matlab.graphics.primitive.Patch')
        set(h.bx,'FaceColor',confs.scattercolors{1});
        set(h.md,'Color',confs.scattercolors{2});
    else
        set(h.bx,'Color',confs.scattercolors{1});
        set(h.md,'Color',confs.scattercolors{1});
    end

    for g = 1:num_groups
       set(h.sc(:,g),'MarkerFaceColor',confs.colors(g,:))
    end
    
    % add a legend based on scatter colors
    if ~isempty(confs.legend)
        h.lg = legend(h.sc(1,:),confs.legend);
    end
else
    
    % add a legend based on box colors
    if ~isempty(confs.legend)
        h.lg = legend(h.bx(1,:),confs.legend);
    end
end

% set ticks and labels
set(gca,'XTick',cpos,'XTickLabels',cpos,'box','off');
if ~isempty(confs.xtlabels)
    set(gca,'XTickLabels',confs.xtlabels,'XTick',cpos);
end  

xlim([gpos(1)-bar_width, gpos(end)+bar_width]); % adjust x-axis margins

end