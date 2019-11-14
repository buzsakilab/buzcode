function [rankStats] = bz_getEventRank(varargin)
% [rankStats] = getEventRank()
%  Get rank order of spikes inside events from previously calculated 
% bz_getEventSpikes structure
%
% INPUTS
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'basepath'	    full path where session is located (default pwd)
%                       e.g. /mnt/Data/buddy140_060813_reo/buddy140_060813_reo
%     'spkEventTimes'	Field to be used is:    
%                   		.EventRel: 2xN cell matrix. In the first row, relative 
%                   		times of spikes for that particular event across all units.
%                   		In the second row, the UID associated to the above spike
%     'templateType'    String that can be:
%                       - 'Pairwise': (default) it computes the rank correlation
%                           of each event pairwise with each other event.
%                       - 'MeanRank': it computes the rank correlation 
%                           of each event against the mean rank
%                           of the rest of the events.
%                       - 'Peak': searchs for the bz_findPlaceFieldsTemplate
%                           output, the X.placeFieldTemplate.mat, loads it
%                           and takes the 'Peak' field, a (# units) x 3 x (# conditions)
%                           matrix, which has bins corresponding to the firing 
%                           map peak for each unit (NaN for units without place
%                           field); second column has the unit ID; third column
%                           has the position of  the unit in the UID vector. The
%                           third dimension contains this similar matrix for the
%                           other conditions.
%                       - 'CenterofMass': searchs for the bz_findPlaceFieldsTemplate
%                           output, the X.placeFieldTemplate.mat, loads it
%                           and takes the 'CenterofMass' field, a (# units) x 
%                           3 x (# conditions) matrix, which has bins corresponding 
%                           to the firing map peak for each unit (NaN for units  
%                           without place field); second column has the unit ID; 
%                           third column has the position of  the unit in the UID 
%                           vector. The third dimension contains this similar matrix 
%                           for the other conditions.
%     'timeSpike'       A string, to determine what time reference of spike:
%                       1. 'first': (default) takes into account the first time
%                           the unit fires 
%                       2. 'mean': takes the mean of the spike time
%     'minUnits'        Minimum number of units parcitipating in each event. 
%     'normalize'       Work with normalized ranks, logical (default: true)
%                       So instead of ranking 1, 4, 5, 2, 3, it will be 0,
%                       0.8, 1.0, 0.2, 0.4.
%     'numRep'          Number of permutation test repetitions 
%                       to compute statistical significance. In each
%                       permutation test, rank of firing units on each
%                       event is randomly permuted (default: 1000)
%     'pvalTest'        pval value under which to test the total statistical 
%                       significance of sequence consistency (default: 0.05)
%     'minCorr'         minimum correlation value to search for rank
%                       clusters, clusters within which rank is highly 
%                       correlated (default: 0.8)
%     'doPlot'	        Make plot, logical (default: true) 
%     'saveMat'	        Saves file, logical (default: true) 
%
%    =========================================================================
%
%
% INPUTS
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'rankStats'	    structure with the following subfields:
%         .corrMean          (double) Mean of rank correlations of all events
%         .corrStd           (double) Standard deviation of rank correlations
%         .corrEvents        (1xN matrix) Rank correlations of each event.
%                            N: number of events
%         .corrMeanShuff     (numRep x 1 matrix) Mean of rank correlations 
%                            of each event for all permutation tests
%         .corrEventsShuff   (numRep x N matrix) Rank correlations 
%                            of each event for all permutation tests
%         .pvalEvents        pvalue per event
%         .pvalTotal         binomial test over 'pvalEvents'
%                            
%                            
%
%  Andrea Navas-Olive, 2019

% TODO:
%  - Have an option to test two group events 
%  - Test With an external template

%% Parse inputs 
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'spkEventTimes',{},@isstruct);
addParameter(p,'templateType','MeanRank',@isstr);
addParameter(p,'timeSpike','first',@isstr);
addParameter(p,'minUnits',5,@isnumeric);
addParameter(p,'normalize',true,@islogical);
addParameter(p,'numRep',1000,@isnumeric);
addParameter(p,'pvalTest',0.05,@isnumeric);
addParameter(p,'minCorr',0.05,@isnumeric);
addParameter(p,'doPlot', true, @islogical);
addParameter(p,'saveMat', true, @islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
spkEventTimes = p.Results.spkEventTimes;
templateType = p.Results.templateType;
normalize = p.Results.normalize;
timeSpike = p.Results.timeSpike;
minUnits = p.Results.minUnits;
numRep = p.Results.numRep;
pvalTest = p.Results.pvalTest;
minCorr = p.Results.minCorr;
doPlot = p.Results.doPlot;
saveMat = p.Results.saveMat;

% Get session info
basename = bz_BasenameFromBasepath(basepath);
load([basepath filesep basename '.sessionInfo.mat']);
% Load default spkEventTimes
if isempty(spkEventTimes)
    spkEventTimes = load([basepath filesep basename '.spkEventTimes.mat']);
    spkEventTimes = spkEventTimes.spkEventTimes;
end
% Relative times of spikes for each particular event across all units.
evtTimes = spkEventTimes.EventRel;

% If template type is 'Peak' or 'CenterofMass' load 'placeFieldTemplate'
if strcmp(templateType,'Peak') || strcmp(templateType,'CenterofMass')
    % External template
    templateExt = load([basepath filesep basename '.placeFieldTemplate.mat']);
    templateExt = templateExt.placeFieldTemplate;
    templateExt = templateExt.(templateType);
else
    templateExt = [];
end


%% Rank matrix
% Create (#events) x (#units) matrix with position of spikes in event. It
% considers just the first spikes of each unit
rankUnits = nan*ones(size(spkEventTimes.UnitEventRel));
for event = 1:length(evtTimes)
    % Take into account just first spike
    if strcmp(timeSpike,'first')
        units = unique(evtTimes{event}(3,:),'stable');
    elseif strcmp(timeSpike,'mean')
        means = [];
        for jj = unique(evtTimes{event}(3,:))
            means  = [means [mean(evtTimes{event}(1,evtTimes{event}(3,:)==jj)); jj]];
        end
        if ~isempty(means)
            means = sortrows(means')';
            units = means(2,:);
        else
            units = [];
        end
    else
        warning('The variable "timeSpike" is invalid');
    end
    nUnits = length(units);
    % Set event as nan if it has no enough units
    if nUnits < minUnits
        rankUnits(:,event) = nan;
    % Rank units
    else
        rankUnits(units,event) = 1:nUnits;
        % If normalize, set ranks between 0 and 1        
        if normalize
            rankUnits(units,event) = rankUnits(units,event) / nUnits;
        end
    end
end

%% Compute correlation

% Compute mean and standard deviation of rank correlation, and correlation
% event by event
[corrMean, corrStd, corrEvents] = compute_rank_correlation(templateType, rankUnits, evtTimes, templateExt);

% Permutation test: random shuffling the ranks of the units that 
% participated in each event and carrying out the entire procedure, 
% culminating in the computation of the mean rank correlation. This is 
% repeated 'numRep' times and the observed mean rank correlation is 
% compared with the distribution of the randomly permuted mean rank 
% correlations
corrMeanShuff = zeros(numRep,1);
corrEventsShuff = zeros(numRep,length(evtTimes));
parfor irep = 1:numRep
    % Suffle each column of rankUnits 
    rankUnitsSuffled = shuffle_by_column(rankUnits,normalize);
    % Rank correlation with suffled rankUnits
    [tmpcorr, ~, tmpcorrTemplate] = compute_rank_correlation(templateType, rankUnitsSuffled, evtTimes, templateExt);
    corrMeanShuff(irep) = tmpcorr;
    corrEventsShuff(irep,:) = tmpcorrTemplate;    
end

% For statistical significance of sequence consistency, get percentage of
% times that the correlation of a certain event has been higher that the
% correlation of that same event in the permutation test
pvalEvents = zeros(1,length(corrEvents));
parfor event = 1:length(corrEvents)
    pvalEvents(event) = sum( corrEvents(event) > corrEventsShuff(:,event) ) / numRep;    
end
% Transform it in a p-value
pvalEvents = 1 - pvalEvents;
% Compute a general p-value, measuring how probable is to get such
% pvalEvents distribution
pvalTotal = binom_test(pvalEvents,pvalTest);


%% Find different sequences

% Events with statistically significant rank sequences
statEvents = find(pvalEvents<=0.05); 

if length(statEvents)>2
    % Perform a pairwise correlation
    corrMat = corr(rankUnits(:,statEvents), 'rows', 'pairwise');
    % Clean correlation matrix (remove nans)
    corrMat(isnan(corrMat)) = 0;
    % Take intro account correlation above 'minCorr'
    corrMat(corrMat <= minCorr) = 0;

    % Perform an agglomerative hierarchical cluster tree
    tree = linkage(corrMat,'complete','correlation');
    % Compute a cutoff (to select the hierarchical level to perform the clusterization)
    treeNonan = tree(~isnan(tree(:,3)),:);
    cutoff = median([treeNonan(end-2,3) treeNonan(end-1,3)]); 
    % Cluster
    groups = cluster(tree,'cutoff',cutoff,'criterion','distance');
    
    % Save sequences
    rankClusters = zeros(size(pvalEvents));
    rankClusters(statEvents) = groups;
else
    rankClusters = zeros(size(pvalEvents));
end



%% Save and plot

% Save statistics
if saveMat
    rankStats.corrMean = corrMean;
    rankStats.corrStd = corrStd;
    rankStats.corrEvents = corrEvents;
    rankStats.corrMeanShuff = corrMeanShuff;
    rankStats.corrEventsShuff = corrEventsShuff;
    rankStats.pvalEvents = pvalEvents;
    rankStats.pvalTotal = pvalTotal;
    rankStats.rankUnits = rankUnits;
    rankStats.rankClusters = rankClusters;
    save([basepath filesep basename '.rankStats.mat'],'rankStats'); 
end
% Plot statistics
if doPlot
    
    % Plot results of permutation test
    figure
    hold on
    [yhist,xhist] = hist(corrEventsShuff,[-1:0.05:1]);
    fill([xhist; flip(xhist)],[mean(yhist,2)+std(yhist')';mean(yhist,2)-std(yhist')']/max(mean(yhist,2)),1,'facecolor',[.7,.9,1],'linestyle','none','DisplayName','')
    plot(xhist,mean(yhist,2)/max(mean(yhist,2)),'b','linewidth',2,'DisplayName','Shuffled data')
    [yhist,xhist] = hist(corrEvents,[-1:0.05:1]);
    plot(xhist,yhist/max(yhist),'k','linewidth',2,'DisplayName','Real data')
    legend()
    xlabel('Correlation values')
    ylabel('Distribution') 
    
    if exist('tree')
        % Plot hierarchical tree 
        figure(),
        dendrogram(tree,'ColorThreshold',cutoff)
        xlabel('Events')

        % Plot rank clusters
        figure('pos',[100,300,1000,400]), hold on
        for group = unique(groups)'
            idxsGroup = groups==group;
            cycsRank = rankUnits(:,idxsGroup);
            meanUnitPos = [];
            for unit = find(sum(isnan(cycsRank),2) <= 0.5*sum(idxsGroup))'
                meanUnitPos = [meanUnitPos [unit; nanmedian(cycsRank(unit,:))]];
            end
            if ~isempty(meanUnitPos)
                meanUnitPos = sortrows(meanUnitPos',2)';
                % Plot in order
                subplot(1,max(groups),group), hold on
                k = -1;
                for unit = meanUnitPos(1,:)
                    color = 'k';
                    plot(nanmedian(rankUnits(unit,idxsGroup)), k,'.k','markersize',30)
                    errorbar(nanmedian(rankUnits(unit,idxsGroup)),k,0,0,nanstd(rankUnits(unit,idxsGroup)),nanstd(rankUnits(unit,idxsGroup)),'color',color)
                    plot(rankUnits(unit,idxsGroup), k + 0.1*randn(1,sum(idxsGroup)),'.')
                    k = k-1;
                end
                set(gca,'ytick',[k+1:-1],'yticklabel',meanUnitPos(1,end:-1:1));
                xlabel('Rank position')
                ylabel('Units')
            end
        end
    end

end


end




% RANK CORRELATION...
function [corrMean, corrStd, corrEvents] = compute_rank_correlation(templateType, spkPos, evtTimes, templateExt)


    % ... WITHOUT EXTERNAL TEMPLATE
    % Template method: a template for each specific ripple or theta event is
    % constructed based on the averaged rank of all units over all other events
    % (i.e., excluding that specific event). Then the rank correlation between
    % each specific event and its template was computed, and averaged over all
    % events
    if strcmp(templateType,'MeanRank')
        corrEvents = zeros(size(evtTimes));
        parfor event = 1:length(evtTimes)
            % Order of units in this event
            rankThisEvent = spkPos(:,event);
            % Mean order of units alon  g the rest of the events (template)
            rankTemplate = nanmean(spkPos(:,[1:event-1,event+1:end]),2);
            % Correlation between these previous variables
            corrEvents(event) = corr(rankThisEvent, rankTemplate, 'rows', 'complete');
        end
        % Mean and standard deviation of rank correlation
        corrMean = nanmean(corrEvents);
        corrStd = nanstd(corrEvents);
    % Pairwise method:
    elseif strcmp(templateType,'Pairwise')
        corrEvents = corr(spkPos, 'rows', 'pairwise');
        corrEvents(find(eye(size(corrEvents)))) = nan;
        corrEvents = nanmean(corrEvents);
        % Mean and standard deviation of rank correlation
        corrMean = nanmean(corrEvents);
        corrStd = nanstd(corrEvents);
    end

    
    % ... WITH EXTERNAL TEMPLATE
    % Template method: external template
    if ismember(templateType,{'Peak','CenterofMass'})
        % Initialize
        corrEvents = zeros(size(evtTimes,2),length(templateExt));
        corrMean = zeros(1,length(templateExt));
        corrStd = zeros(1,length(templateExt));
        for iCond = 1:length(templateExt)
            parfor event = 1:size(evtTimes,2)
                % Order of units in this event
                rankThisEvent = spkPos(:,event);
                % External template (we need to transform it to an array
                % similar to rankThisEvent: a (# units) x 1 array, where in
                % the position of each unit the rank is written 
                rankTemplate = nan*ones(size(rankThisEvent));
                rankTemplate(templateExt{iCond}(:,3)) = [1:length(templateExt{iCond}(:,3))];
                % Normalize
                if normalize
                    rankTemplate = rankTemplate / length(templateExt{iCond}(:,3));
                end
                % Correlation between these previous variables
                corrEvents(event,iCond) = corr(rankThisEvent, rankTemplate, 'rows', 'complete');
            end
        end
        % Mean and standard deviation of rank correlation
        corrMean(iCond) = nanmean(corrEvents);
        corrStd(iCond) = nanstd(corrEvents);
    end
end

% SUFFLE RANK IN EACH EVENT
function A = shuffle_by_column(A,normalize)
    % Nan temporal vector
    tmpNan = nan * ones(size(A,1),1);
    % For each event...
    parfor ii = 1:size(A,2)
        % Find units participating
        rankAi = find(~isnan(A(:,ii)));
        % Shuffle ranks
        tmp = tmpNan; 
        tmp(rankAi) = randperm(length(rankAi));
        if normalize
            tmp = tmp / length(rankAi);
        end
        % Update matrix
        A(:,ii) = tmp;
    end
end

% BINOMIAL TEST
function p = binom_test(ps,alpha)
%   GLOBAL_P=BINOM_TEST(PS,ALPHA) where PS are the p-values of independent tests
%   and ALPHA is 0.05 per default.
    if nargin<2 || isempty(alpha)
        alpha=.05;
    end
    ps(isnan(ps)) = [];
    % Number of significant pvals
    k = sum(ps<=alpha);
    % Total number of pvals
    n = numel(ps);
    % Sum of: (alpha^q) * (1-alpha)^(n-q) * (n!/(q! (n-q)!)), q from k to n 
    p = nansum( feval( @(q) arrayfun( @(k) nchoosek(n,k), q) * ( alpha.^q.*(1-alpha).^(n-q) )', k:n ));
end
