function [counts,indices] = RadialMaze(data,baitedArms,varargin)

%RadialMaze - Compute simple statistics for the 8-arm radial maze task.
%
% Compute the number of reference and working memory errors, rewards found, etc.
%
%  USAGE
%
%    [counts,indices] = RadialMaze(data,baitedArms,<options>)
%
%    data           list of [rat day trial configuration group arm] tuples
%    baitedArms     list of baited arms (from 1 to 8) for each (rat,configuration)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'plots'       analyses to plot (default = {'individuals','groups'})
%     'measures'    behavioral measures to compute/plot (default = {'perf'})
%     'anova'       perform ANOVAs for selected measures (default = 'on')
%    =========================================================================
%
%  OUTPUT
%
%
%    counts 	a N x 12 matrix, where each line lists one case
%                (rat,day,trial,config,group,measure1,...,measureN)
%    				where the measures are nRefErrors, nWorkErrors, nRewardsFound,
%              nVisits, performance, nPenaltyErrors, nRefErrors, and nAdjustedErrors
%    indices 	in the same format (N x 12 matrix) as counts
%
%
%  NOTE
%
%    The option 'plots' determines whether individual and/or group analyses should
%    be plotted. The option 'measures' determines which kinds of analyses should be
%    plotted: possible values include 'ref','work','rewards','visits','perf','penalty',
%    and 'adjusted'.

% Copyright (C) 2008-2011 by Gabrielle Girardeau & Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% We use a few globals to simplify the code
global showANOVAs;
global measures plots;
global nTrials nDays nArms nConfigurations nRats nGroups;
global nVisitsDayOne;
global p0 pInf;
global RAT DAY TRIAL CONF GROUP ARM REF WORK REWARDS VISITS PERF PENALTY ADJUSTED;
global alpha;

alpha = 0.05; % Significance level for complementary tests

% Column names for the 'data', 'counts' and 'indices' matrices
RAT = 1;
DAY = 2;
TRIAL = 3;
CONF = 4;
GROUP = 5;

% Column names for the 'data' matrix
ARM = 6;

% Column names for the 'counts' and 'indices' matrices
REF = 6;
WORK = 7;
REWARDS = 8;
VISITS = 9;
PERF = 10;
PENALTY = 11;
ADJUSTED = 12;

% Constants
nDays = max(data(:,DAY));
nTrials = 3;
nArms = 8;
nGroups = max(data(:,GROUP));
nConfigurations = max(data(:,CONF));
nRats = length(unique(data(:,RAT)));

% Used later to estimate chance levels for performance
nVisitsDayOne = round(Accumulate(data(:,DAY),1)/nRats/nTrials);
nVisitsDayOne = nVisitsDayOne(1);

% Default values
allMeasures = {'ref','work','rewards','visits','perf','penalty','adjusted'};
allPlots = {'groups','individuals'};
defaultMeasures = {'perf'};
defaultPlots = {'groups','individuals'};
listMeasures = false;
listPlots = false;
measures = {};
plots = {};
showANOVAs = true;

if nargin < 2,
	error('Incorrect number of parameters (type ''help <a href="matlab:help RadialMaze">RadialMaze</a>'' for details).');
end

if mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help RadialMaze">RadialMaze</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+firstIndex) ' is not a property (type ''help <a href="matlab:help RadialMaze">RadialMaze</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'measures',
			p = lower(varargin{i+1});
			if ~isa(p,'cell'),
				error('Value associated with property ''measures'' is not a cell array of strings (type ''help <a href="matlab:help RadialMaze">RadialMaze</a>'' for details).');
			end
			for j = 1:length(p),
				if ~isstring_FMAT(p{j},allMeasures{:}),
					error('Incorrect value for property ''measures'' (type ''help <a href="matlab:help RadialMaze">RadialMaze</a>'' for details).');
				end
			end
			measures = {measures{:},p{:}};
			listMeasures = true;
		case 'plots',
			f = lower(varargin{i+1});
			if ~isa(f,'cell'),
				error('Value associated with property ''plots'' is not a cell array of strings (type ''help <a href="matlab:help RadialMaze">RadialMaze</a>'' for details).');
			end
			for j = 1:length(f),
			if ~isstring_FMAT(f{j},allPlots{:}),
					error('Incorrect value for property ''plots'' (type ''help <a href="matlab:help RadialMaze">RadialMaze</a>'' for details).');
				end
			end
			plots = {plots{:},f{:}};
			listPlots = true;
		case 'anova',
			doIt = lower(varargin{i+1});
			if ~isstring_FMAT(doIt,'on','off'),
				error('Incorrect value for property ''anova'' (type ''help <a href="matlab:help RadialMaze">RadialMaze</a>'' for details).');
			end
			showANOVAs = strcmp(doIt,'on');
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help RadialMaze">RadialMaze</a>'' for details).']);
	end
end

if ~listMeasures, measures = defaultMeasures; end
if ~listPlots, plots = defaultPlots; end

% Compute stats for each rat and each configuration
counts = Counts(data,baitedArms);
indices = ComputeIndices(counts);

EstimatePerformanceChanceLevels;
%  p0 = 0.481869891541284; pInf = 0.306735773403144;

h = ComputeANOVAs(indices);

if any(ismember('individuals',plots)),
	PlotIndividualStats(counts);
end
if any(ismember('groups',plots)),
	PlotGroupStats(indices,h);
end

%-----------------------------------------------------------------------------------------------------------------
%  Count errors (etc.)
%
%  INPUT
%    data, baitedArms
%
%  OUTPUT
%    counts, a N x 12 matrix, where each line lists one case (rat,day,trial,config,group,measure1,...,measureN)
%    where the measures are nRefErrors, nWorkErrors, nRewardsFound, nVisits, performance, nPenaltyErrors,
%    nRefErrors, and nAdjustedErrors
%-----------------------------------------------------------------------------------------------------------------

function counts = Counts(data,baitedArms)

global nTrials nDays nArms nConfigurations nGroups;
global RAT DAY TRIAL CONF GROUP ARM REF WORK REWARDS VISITS PERF PENALTY ADJUSTED;

nRats = max(data(:,RAT));

% Create output array and prefill it with NaNs so missing data will be readily visible
% 5 variables (rat, day, trial, config, group) and 7 measures (ref, work, etc.)
%  counts = nan(nRats*nDays*nTrials*nConfigurations,5+7);

% Let's go...
line = 1;
for rat = 1:nRats,

	% Which group does this rat belong to?
	group = data(data(:,RAT)==rat,GROUP);
	group = group(1);

	for configuration = 1:nConfigurations,

		% Determine baited/unbaited arms for this rat and this configuration
		b = baitedArms(rat,configuration,:);
		baitedArm1 = b(1);
		baitedArm2 = b(2);
		baitedArm3 = b(3);
		unbaitedArms = setdiff(1:8,[baitedArm1 baitedArm2 baitedArm3]);
		unbaitedArm1 = unbaitedArms(1);
		unbaitedArm2 = unbaitedArms(2);
		unbaitedArm3 = unbaitedArms(3);
		unbaitedArm4 = unbaitedArms(4);
		unbaitedArm5 = unbaitedArms(5);

		for day = 1:nDays,
			for trial = 1:nTrials,

				% Select the block of data corresponding the current rat, config, trial and day
				block = data(data(:,RAT)==rat&data(:,CONF)==configuration&data(:,DAY)==day&data(:,TRIAL)==trial,:);

				% Total number of visits
				nVisits = size(block,1);

				% Missing data: skip
				if nVisits == 0,
					continue;
				end

				% Set values of the variables for this line in 'counts'
				counts(line,1:REF-1) = [rat,day,trial,configuration,group];

				% Initially set all measures to zero (except number of visits!)
				counts(line,REF:ADJUSTED) = zeros(1,7);
				counts(line,VISITS) = nVisits;

				% No arm visited yet
				visitedBaitedArm1 = false;
				visitedBaitedArm2 = false;
				visitedBaitedArm3 = false;
				visitedUnbaitedArm1 = false;
				visitedUnbaitedArm2 = false;
				visitedUnbaitedArm3 = false;
				visitedUnbaitedArm4 = false;
				visitedUnbaitedArm5 = false;

				for visit = 1:nVisits,
					% Test successively visited arms to count the number of working memory errors, reference memory errors, etc.
					% In particular, check for repeated visits in baited arms (returning in an already visited baited arm is considered a working memory error)
					currentArm = block(visit,ARM);
					switch currentArm,
						case baitedArm1,
							if visitedBaitedArm1,
								% Returning in an already visited baited arm is considered as a working memory error
								counts(line,WORK) = counts(line,WORK)+1;
							else
								visitedBaitedArm1 = true;
								counts(line,REWARDS) = counts(line,REWARDS)+1;
							end
						case baitedArm2,
							if visitedBaitedArm2,
								counts(line,WORK) = counts(line,WORK)+1;
							else
								visitedBaitedArm2 = true;
								counts(line,REWARDS) = counts(line,REWARDS)+1;
							end
						case baitedArm3,
							if visitedBaitedArm3,
								counts(line,WORK) = counts(line,WORK)+1;
							else
								visitedBaitedArm3 = true;
								counts(line,REWARDS) = counts(line,REWARDS)+1;
							end
						case unbaitedArm1,
							if visitedUnbaitedArm1,
								counts(line,WORK) = counts(line,WORK)+1;
							else
								visitedUnbaitedArm1 = true;
								counts(line,REF) = counts(line,REF)+1;
							end
						case unbaitedArm2,
							if visitedUnbaitedArm2,
								counts(line,WORK) = counts(line,WORK)+1;
							else
								visitedUnbaitedArm2 = true;
								counts(line,REF) = counts(line,REF)+1;
							end
						case unbaitedArm3,
							if visitedUnbaitedArm3,
								counts(line,WORK) = counts(line,WORK)+1;
							else
								visitedUnbaitedArm3 = true;
								counts(line,REF) = counts(line,REF)+1;
							end
						case unbaitedArm4,
							if visitedUnbaitedArm4,
								counts(line,WORK) = counts(line,WORK)+1;
							else
								visitedUnbaitedArm4 = true;
								counts(line,REF) = counts(line,REF)+1;
							end
						case unbaitedArm5,
							if visitedUnbaitedArm5,
								counts(line,WORK) = counts(line,WORK)+1;
							else
								visitedUnbaitedArm5 = true;
								counts(line,REF) = counts(line,REF)+1;
							end
					end
				end
				counts(line,PERF) = counts(line,REWARDS)/counts(line,VISITS);
				% Adjusted errors = visits in non-rewarded arms + unvisited baited arms
				counts(line,ADJUSTED) = counts(line,REF)+(3-counts(line,REWARDS));
				% Penalty errors = 5 if the rat did not find all the rewards
				if counts(line,REWARDS)~= 3,
					counts(line,PENALTY) = 5;
				else
					counts(line,PENALTY) = counts(line,REF);
				end
				line = line + 1;
			end
		end
	end
end


%-----------------------------------------------------------------------------------------------------------------
%  Transform counts/proportions for ANOVA (square root transform, arcsine square root transform)
%
%  INPUT
%    counts
%
%  OUTPUT
%    indices, in the same format (N x 12 matrix) as counts
%-----------------------------------------------------------------------------------------------------------------

function indices = ComputeIndices(counts)

global RAT DAY TRIAL CONF GROUP ARM REF WORK REWARDS VISITS PERF PENALTY ADJUSTED;

indices = counts;
indices(:,[REF WORK REWARDS VISITS PENALTY ADJUSTED]) = ComputeCountIndex(counts(:,[REF WORK REWARDS VISITS PENALTY ADJUSTED]));
indices(:,PERF) = ComputeProportionIndex(counts(:,PERF));

function index = ComputeCountIndex(count)

index = sqrt(count+1);

function index = ComputeProportionIndex(proportion)

index = asin(sqrt(proportion))/(pi/2);

%-----------------------------------------------------------------------------------------------------------------
%  Individual plots: plot error counts for each rat in each configuration
%
%  INPUT
%    counts
%-----------------------------------------------------------------------------------------------------------------

function PlotIndividualStats(counts)

global nTrials nDays nConfigurations measures plots testRats controlRats;
global RAT DAY TRIAL CONF GROUP ARM REF WORK REWARDS VISITS PERF PENALTY ADJUSTED;

if any(ismember('ref',measures)),
	DoPlotIndividualStats(counts,REF,'# ref mem errors',[0 5]);
end
if any(ismember('work',measures)),
	DoPlotIndividualStats(counts,WORK,'# work mem errors',[0 5]);
end
if any(ismember('visits',measures)),
	DoPlotIndividualStats(counts,VISITS,'# visits',[0 10]);
end
if any(ismember('rewards',measures)),
	DoPlotIndividualStats(counts,REWARDS,'# rewards',[0 3]);
end
if any(ismember('perf',measures)),
	DoPlotIndividualStats(counts,PERF,'performance',[0 1]);
end
if any(ismember('penalty',measures)),
	DoPlotIndividualStats(counts,PENALTY,'# penalty errors',[0 5]);
end
if any(ismember('adjusted',measures)),
	DoPlotIndividualStats(counts,ADJUSTED,'# adjusted errors',[0 5]);
end


function DoPlotIndividualStats(counts,measure,measureStr,yLim)

global nRats nConfigurations nDays;
global RAT DAY TRIAL CONF GROUP ARM REF WORK REWARDS VISITS PERF PENALTY ADJUSTED;

colors = 'rbkgcmy'; % Hopefully there will not be more than 7 groups!
for configuration = 1:nConfigurations,
	figure;
	for rat = 1:nRats,
		SquareSubplot(nRats,rat);
		d = counts(counts(:,RAT)==rat&counts(:,CONF)==configuration,[DAY measure GROUP]);
		if isempty(d), continue; end
		y = Accumulate(d(:,1),d(:,2));
		n = Accumulate(d(:,1),1);
		plot(y./n,colors(d(1,3)));
		xlabel(['Rat #' int2str(rat) ' - Configuration #' int2str(configuration)]);
		ylabel(measureStr);
		xlim([0 nDays+1]);set(gca,'xtick',1:nDays);ylim(yLim);
	end
end

%-----------------------------------------------------------------------------------------------------------------
%  Group plots: plot error indices for all rats in each configuration
%
%  INPUT
%    indices
%    h   post-hoc ttest results from ANOVAs
%-----------------------------------------------------------------------------------------------------------------

function PlotGroupStats(indices,h)

global measures;
global RAT DAY TRIAL CONF GROUP ARM REF WORK REWARDS VISITS PERF PENALTY ADJUSTED;

% Reference memory errors
if any(ismember('ref',measures)),
	DoPlotGroupStats(indices,REF,'index ref mem errors',[0 5],h);
end
if any(ismember('work',measures)),
	DoPlotGroupStats(indices,WORK,'index work mem errors',[0.5 2],h);
end
if any(ismember('rewards',measures)),
	DoPlotGroupStats(indices,REWARDS,'index rewards errors',[0.3 1],h);
end
if any(ismember('visits',measures)),
	DoPlotGroupStats(indices,VISITS,'index visits errors',[0.3 1],h);
end
if any(ismember('perf',measures)),
	DoPlotGroupStats(indices,PERF,'index performance errors',[0.3 1],h);
end
if any(ismember('penalty',measures)),
	DoPlotGroupStats(indices,PENALTY,'index penalty errors',[0.3 1],h);
end
if any(ismember('adjusted',measures)),
	DoPlotGroupStats(indices,ADJUSTED,'index adjusted errors',[0.3 1],h);
end


function DoPlotGroupStats(indices,measure,measureStr,yLim,h)

global nTrials nDays nConfigurations nRats nGroups;
global RAT DAY TRIAL CONF GROUP ARM REF WORK REWARDS VISITS PERF PENALTY ADJUSTED;
global p0 pInf;

colors = 'rbkgcmy'; % Hopefully there will not be more than 7 groups!
for configuration = 1:nConfigurations,
	figure;hold on;
	for group = 1:nGroups,
		d = indices(indices(:,GROUP)==group&indices(:,CONF)==configuration,[DAY measure]);
		if isempty(d), continue; end
		[m,s] = Accumulate(d(:,1),d(:,2));
		n = Accumulate(d(:,1),1);
		sem = s./sqrt(n);
		p = errorbar(m,sem,colors(group));
		set(p,'marker','.','markersize',24);
		xlabel(['Configuration #' int2str(configuration)]);
		ylabel(['index ' measureStr]);
		xlim([0 nDays+1]);set(gca,'xtick',1:nDays);ylim(yLim);

		if ~isempty(h),
			x = find(h(:,1+group*2,measure));
			plot(x,(m(x)-sem(x))-.025,'k+','markersize',10);
			if group == 1 & configuration == 1,
				x = find(h(:,1,measure));
				plot(x,(m(x)+sem(x))+.025,'k*','markersize',10);
			end
		end
	end
	if measure == PERF,
		PlotIntervals([pInf,p0],'rectangles','h');
	end
end


%-----------------------------------------------------------------------------------------------------------------
%  Compute ANOVAs
%
%  INPUT
%    indices
%-----------------------------------------------------------------------------------------------------------------

function h = ComputeANOVAs(indices)

global measures;
global RAT DAY TRIAL CONF GROUP ARM REF WORK REWARDS VISITS PERF PENALTY ADJUSTED;

if any(ismember('ref',measures)),
	h(:,:,REF) = DoComputeANOVAs(indices,REF,'ref mem errors');
end
if any(ismember('work',measures)),
	h(:,:,WORK) = DoComputeANOVAs(indices,WORK,'work mem errors');
end
if any(ismember('visits',measures)),
	h(:,:,VISITS) = DoComputeANOVAs(indices,VISITS,'visits');
end
if any(ismember('rewards',measures)),
	h(:,:,REWARDS) = DoComputeANOVAs(indices,REWARDS,'rewards');
end
if any(ismember('perf',measures)),
	h(:,:,PERF) = DoComputeANOVAs(indices,PERF,'performance');
end
if any(ismember('penalty',measures)),
	h(:,:,PENALTY) = DoComputeANOVAs(indices,PENALTY,'penalty errors');
end
if any(ismember('adjusted',measures)),
	h(:,:,ADJUSTED) = DoComputeANOVAs(indices,ADJUSTED,'adjusted errors');
end


function h = DoComputeANOVAs(indices,measure,measureStr)

global showANOVAs;
global RAT DAY TRIAL CONF GROUP ARM REF WORK REWARDS VISITS PERF PENALTY ADJUSTED;

[p,t,stats,terms,h] = ComputeTwoWayANOVA(indices,measure);
if showANOVAs,
	f = statdisptable(t,['ANOVA for ' measureStr],'Analysis of Variance');
	set(f,'numbertitle','off');
end

%-----------------------------------------------------------------------------------------------------------------
%  Compute day x group ANOVA on indices for a given measure (ref errors, work errors, etc.)
%  !!! Only for configuration 1 !!!
%
%  INPUT
%    indices
%
%  OUTPUT
%    p,t,stats,terms   ANOVA results
%    h                 = [H0,p0,H1,p1,...Hn,pn] where [H0,p0] is the result of a t-test comparing groups 1 and 2
%                      (one row per day), and [Hi,pi] is the result of a t-test comparing group i to the upper
%                      chance level
%-----------------------------------------------------------------------------------------------------------------

function [p,t,stats,terms,h] = ComputeTwoWayANOVA(indices,measure)

global RAT DAY TRIAL CONF GROUP ARM REF WORK REWARDS VISITS PERF PENALTY ADJUSTED;
global p0 pInf;
global nDays nGroups alpha;

i = indices(:,CONF) == 1;
[p,t,stats,terms] = anovan(indices(i,measure),{indices(i,GROUP) indices(i,DAY)},'model','full','varnames',{'group','day'},'display','off');

h = zeros(nDays,2+2*nGroups);
switch measure,
	case { REF, WORK, VISITS, PENALTY, ADJUSTED },
		tail = 'right';
	case { PERF, REWARDS },
		tail = 'left';
end

% Compare groups 1 and 2, and each group to chance level
if all(p<0.05),
	for day = 1:nDays,
		d1 = indices(i&indices(:,DAY)==day&indices(:,GROUP)==1,measure);
		d2 = indices(i&indices(:,DAY)==day&indices(:,GROUP)==2,measure);
		if ~isempty(d1) && ~isempty(d2),
%  			keyboard
%  			[h(day,1),h(day,2)] = ttest2(d1,d2,alpha,tail,'equal');
			[h(day,1),h(day,2)] = ttest2(d1,d2,alpha,tail);
		end
		for group = 1:nGroups,
			d0 = indices(i&indices(:,DAY)==day&indices(:,GROUP)==group,measure);
			if ~isempty(d0),
				[h(day,1+2*group),h(day,2+2*group)] = ttest(d0,p0,alpha,'right');
			end
		end
	end
end

%-----------------------------------------------------------------------------------------------------------------
%  DEPRECATED CODE
%-----------------------------------------------------------------------------------------------------------------
%  Compute repeated measures ANOVA on indices for a given measure (ref errors, work errors, etc.)
%  This can be used to compare e.g. test rats and control rats on the same configuration, or test rats on both
%  configurations
%
%  INPUT
%    data1           'indices' for group 1
%    data2           'indices' for group 2
%-----------------------------------------------------------------------------------------------------------------

function ComputeRepeatedMeasuresANOVA(data1,data2)

global nTrials nDays;

group1 = ones(size(data1));
n = size(data1,2);
days1 = repmat((1:nDays)',1,n);
trials1 = repmat((1:nTrials),nDays,n/3);
rats1 = repmat(ceil((1:n)/3),nDays,1);

group2 = 2*ones(size(data2));
n = size(data2,2);
days2 = repmat((1:nDays)',1,n);
trials2 = repmat((1:nTrials),nDays,n/3);
rats2 = repmat(ceil((1:n)/3),nDays,1);

data = [data1(:);data2(:)];
groups = [group1(:);group2(:)];
days = [days1(:);days2(:)];
trials = [trials1(:);trials2(:)];
rats = [rats1(:);rats2(:)];

nans = isnan(data(:,1));
data(nans,:) = [];
groups(nans,:) = [];
days(nans,:) = [];
trials(nans,:) = [];
rats(nans,:) = [];

d = [data groups rats days];

RMAOV1MS(d);

%-----------------------------------------------------------------------------------------------------------------
%  Estimate chance levels for performance
%-----------------------------------------------------------------------------------------------------------------

function EstimatePerformanceChanceLevels

global p0 pInf;

p0 = mean(bootstrp(10000,@SimulateOneTrial,0));
pInf = mean(bootstrp(10000,@SimulateOneTrial,Inf));

% nWME = # working memory errors
% If nWME is Inf, then the rat can choose any arm on each trial, independent of previous choices
% Otherwise, the performance is computed in two steps: 1) the rat never makes working memory errors
% and 2) a fixed number of working memory errors are added to the performance computed in 1).

function index = SimulateOneTrial(nWME)

global nVisitsDayOne;

if isinf(nWME),
	% Draw n integers from 1 to 8 (= arm number)
	n = nVisitsDayOne;
	r = floor(Clip(rand(n,1),0,1-eps)*8)+1;
	% How many 'arms' came out from the list [1 2 3]?
	found = any(r==1) + any(r==2) + any(r==3);
else
	% Find arm #1
	r = randperm(8);
	n1 = find(r==1|r==2|r==3);n1 = n1(1);
	% Find arm #2
	r = randperm(8-n1);
	n2 = find(r==1|r==2);n2 = n2(1);
	% Find arm #3
	r = randperm(8-n1-n2);
	n3 = find(r==1);
	found = 3;
	n = n1+n2+n3+nWME;
end
% Corresponding proportion index
index = ComputeProportionIndex(found/n);
