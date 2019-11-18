function [rankStats] = bz_RankOrder(varargin)
% [rankStats] = RankOrder()
%  Get rank order of spikes inside events from previously calculated 
% bz_getSpikesRank structure. The 'corrEvents' output field structure shows
% rank correlation of each event with the selected 'templateType', and the
% 'pvalEvents' field states how significantly different from chance is that
% correlation, so presumably those events with low 'pvalEvents' will have a
% rank order that repeats over time. In the case that there is more than 
% one rank sequence, the output field structure 'rankClusters' will number
% the events equally if they have similar sequences.
%
% INPUTS
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'basepath'	    full path where session is located (default pwd)
%                       E.g. /mnt/Data/buddy140_060813_reo/buddy140_060813_reo
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
%                           and takes the 'Peak' field, a 1 x (# conditions) 
%                           cell array. Within each cell there is a (# units) x 3
%                           matrix, which has bins corresponding to the firing 
%                           map peak for each unit (NaN for units without place
%                           field); second column has the unit ID; third column
%                           has the position of  the unit in the UID vector. The
%                           rest of the cells contain this similar matrix for the
%                           other conditions.
%                       - 'CenterofMass': searchs for the bz_findPlaceFieldsTemplate
%                           output, the X.placeFieldTemplate.mat, loads it
%                           and takes the 'CenterofMass', a 1 x (# conditions) 
%                           cell array. Within each cell there is a (# units) x 3
%                           matrix, which has bins corresponding to the firing 
%                           map peak for each unit (NaN for units without place 
%                           field); second column has the unit ID; third column
%                           has the position of  the unit in the UID vector. 
%                           The rest of the cells contain this similar matrix 
%                           for the other conditions.
%     'eventIDs'        A (#events)-length array of 0s and 1s that indicates
%                       which event belongs to one of two different groups 
%                       of events. By default all events will belong to a
%                       single group, and the rank correlation of each event
%                       will be performed against the rest of the events.
%                       However, if an 'eventIDs' array establishes two
%                       different groups of events (e.g. spontaneous vs induced
%                       ripples), then the rank correlation of events in 
%                       one group will be performed against the rank of the
%                       other group.
%     'timeSpike'       A string, to determine what time reference of spike:
%                       1. 'first': (default) takes into account the first time
%                           the unit fires 
%                       2. 'mean': takes the mean of the spike time
%     'minUnits'        Minimum number of units parcitipating in each event.
%                       Events with less than 'minUnits' are not considered
%                       for correlation.
%     'normalize'       Work with normalized ranks, logical (default: true)
%                       So instead of ranking 1, 4, 5, 2, 3, it will be 0,
%                       0.8, 1.0, 0.2, 0.4.
%     'numRep'          Number of permutation test repetitions 
%                       to compute statistical significance. In each
%                       permutation test, rank of firing units on each
%                       event is randomly permuted (default: 1000)
%     'pvalTest'        pval value under which to test the total statistical 
%                       significance of sequence consistency (default: 0.05)
%     'minCorr'         minimum correlation value to search for different 
%                       rank clusters (that can be interpreted as different
%                       sequences), clusters within which rank is highly 
%                       correlated (default: 0.8)
%     'doPlot'	        Make plots, logical (default: true) 
%     'saveMat'	        Saves file, logical (default: true) 
%
%    =========================================================================
%
%
% OUTPUTS
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'rankStats'	    structure with the following subfields:
%         .corrMean          (#conditions x 1) Mean of rank correlations 
%                            of all events. If there are several conditions
%                            each row is the mean correlation of events
%                            against the template given by each condition.
%         .corrStd           (#conditions x 1) Standard deviation of rank 
%                            of all events. If there are several conditions
%                            each row is the mean correlation of events
%                            against the template given by each condition.
%         .corrEvents        (#conditions x #events matrix) Rank correlation
%                            of each event. Same as with 'corrMean' happens
%                            for each row.
%         .corrMeanShuff     (numRep x #conditions matrix) Mean of rank 
%                            correlations of each event for all permutation 
%                            tests. Same as with 'corrEvents' happens for 
%                            each column. 
%         .corrEventsShuff   (numRep x #events x #conditions matrix) Rank 
%                            correlations of each event for all permutation
%                            tests. Same as with 'corrEvents' happens for 
%                            the third dimension. 
%         .pvalEvents        (#conditions x 1) P-value per event. If there
%                            are several conditions each row is the p-value
%                            of the correlation between each event rank
%                            against the external template given by each
%                            condition.
%         .pvalTotal         binomial test over 'pvalEvents'
%         .rankUnits         (#units x #events matrix) Rank of units for
%                            each event. This is the matrix that has been
%                            used to compute the correlation.
%                            E.g.:  A  [ 0.0  0.4  0.6  nan       0.6 
%                                   B    0.2  0.2  nan  nan       nan
%                                   C    0.4  0.0  nan  nan       nan
%                                   D    0.6  0.6  1.0  nan  ...  1.0
%                                   E    0.8  0.8  0.3  nan       0.3
%                                   F    1.0  1.0  nan  nan       nan ]
%         .rankClusters      (1 x #events array) Different sequences are
%                            searched using an agglomerative hierarchical
%                            cluster tree. This array shows which of the
%                            different found sequences has each event.
%                            E.g.:    [  1    1    2   nan  ...   2  ]
%                                  where 1 is sequence A B C D E F
%                                        2 is sequence E A D
%                                  
%
%  See also bz_getSpikesRank, bz_findPlaceFieldsTemplate
%
%
%
%  Andrea Navas-Olive, 2019


%% Parse inputs 
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'spkEventTimes',{},@isstruct);
addParameter(p,'templateType','MeanRank',@isstr);
addParameter(p,'eventIDs',1,@isnumeric);
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
eventIDs = p.Results.eventIDs;
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
[corrMean, corrStd, corrEvents] = compute_rank_correlation(templateType, rankUnits, evtTimes, templateExt, eventIDs);

% Permutation test: random shuffling the ranks of the units that 
% participated in each event and carrying out the entire procedure, 
% culminating in the computation of the mean rank correlation. This is 
% repeated 'numRep' times and the observed mean rank correlation is 
% compared with the distribution of the randomly permuted mean rank 
% correlations
nCond = size(corrEvents,1);
corrMeanShuff = zeros(numRep,nCond);
corrEventsShuff = zeros(numRep,length(evtTimes),nCond);
parfor irep = 1:numRep
    % Suffle each column of rankUnits 
    rankUnitsSuffled = shuffle_by_column(rankUnits,normalize);
    % Rank correlation with suffled rankUnits
    [tmpcorr, ~, tmpcorrTemplate] = compute_rank_correlation(templateType, rankUnitsSuffled, evtTimes, templateExt, eventIDs);
    corrMeanShuff(irep,:) = tmpcorr';
    corrEventsShuff(irep,:,:) = tmpcorrTemplate';    
end

% For statistical significance of sequence consistency, get percentage of
% times that the correlation of a certain event has been higher that the
% correlation of that same event in the permutation test
pvalEvents = zeros(nCond,length(corrEvents));
parfor event = 1:length(corrEvents)
    pvalEvents(:,event) = sum( corrEvents(:,event) > squeeze(corrEventsShuff(:,event,:))',2 ) / numRep;    
end
% Transform it in a p-value
pvalEvents = 1 - pvalEvents;
% Compute a general p-value, measuring how probable is to get such
% pvalEvents distribution
pvalTotal = binom_test(pvalEvents,pvalTest);


%% Find different sequences

if any(ismember(templateType,{'Pairwise','MeanRank'}))
    % Events with statistically significant rank sequences
    statEvents = find(pvalEvents<=0.05); 

    try 
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
    catch
        rankClusters = zeros(size(pvalEvents));
    end
else
    rankClusters = [];
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
    
    if exist('groups')
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
function [corrMean, corrStd, corrEvents] = compute_rank_correlation(templateType, rankUnits, evtTimes, templateExt, eventIDs)


    % ... WITHOUT EXTERNAL TEMPLATE
    
    % Template method: a template for each specific ripple or theta event is
    %   constructed based on the averaged rank of all units over all other events
    %   (i.e., excluding that specific event). Then the rank correlation between
    %   each specific event and its template was computed, and averaged over all
    %   events
    if strcmp(templateType,'MeanRank')
        corrEvents = zeros(size(evtTimes));
        % - If there are no two groups of events
        if length(unique(eventIDs))==1
            parfor event = 1:length(evtTimes)
                % Order of units in this event
                rankThisEvent = rankUnits(:,event);
                % Mean order of units along the rest of the events (template)
                rankTemplate = nanmean(rankUnits(:,[1:event-1,event+1:end]),2);
                % Correlation between these previous variables
                corrEvents(event) = corr(rankThisEvent, rankTemplate, 'rows', 'complete');
            end
        % - If there are two groups of events, test this event against 
        %   the meank rank of the other group
        else
            parfor event = 1:length(evtTimes)
                % Order of units in this event
                rankThisEvent = rankUnits(:,event);
                % Mean order of units along the rest of the events (template)
                rankTemplate = nanmean(rankUnits(:,eventIDs==(1-eventIDs(event))),2);
                % Correlation between these previous variables
                corrEvents(event) = corr(rankThisEvent, rankTemplate, 'rows', 'complete');
            end
        end
        % Mean and standard deviation of rank correlation
        corrMean = nanmean(corrEvents);
        corrStd = nanstd(corrEvents);
        
    % Pairwise method:
    elseif strcmp(templateType,'Pairwise')
        corrEvents = corr(rankUnits, 'rows', 'pairwise');
        corrEvents(eye(size(corrEvents))==1) = nan;
        % Mean and standard deviation of rank correlation...
        % - If there are no two groups of events
        if length(unique(eventIDs))==1
            corrEvents = nanmean(corrEvents);
        % - If there are two groups of events, test this event against 
        %   the meank rank of the other group
        else
            corrEventsMean = zeros(size(evtTimes));
            parfor event = 1:length(evtTimes)
                corrEventsMean(event) = nanmean(corrEvents(event,eventIDs==(1-eventIDs(event))));
            end
            corrEvents = corrEventsMean;
        end
        % Mean and standard deviation of rank correlation
        corrMean = nanmean(corrEvents);
        corrStd = nanstd(corrEvents);
    end

    
    % ... WITH EXTERNAL TEMPLATE
    % Template method: external template
    if ismember(templateType,{'Peak','CenterofMass'})
        % Initialize
        corrEvents = zeros(length(templateExt),size(evtTimes,2));
        corrMean = zeros(length(templateExt),1);
        corrStd = zeros(length(templateExt),1);
        for iCond = 1:length(templateExt)
            parfor event = 1:size(evtTimes,2)
                % Order of units in this event
                rankThisEvent = rankUnits(:,event);
                % External template (we need to transform it to an array
                % similar to rankThisEvent: a (# units) x 1 array, where in
                % the position of each unit the rank is written 
                rankTemplate = nan*ones(size(rankThisEvent));
                rankTemplate(templateExt{iCond}(:,3)) = templateExt{iCond}(:,1);
                % Normalize
                if max(max(rankUnits))==1
                    rankTemplate = (rankTemplate - nanmin(templateExt{iCond}(:,1)))/ (nanmax(templateExt{iCond}(:,1))- nanmin(templateExt{iCond}(:,1)));
                end
                % Correlation between these previous variables
                corrEvents(iCond,event) = corr(rankThisEvent, rankTemplate, 'rows', 'complete');
            end
        end
        % Mean and standard deviation of rank correlation
        corrMean = nanmean(corrEvents,2);
        corrStd = nanstd(corrEvents,[],2);
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
