function turnDistribution = RadialMazeTurns(data,mode)

%RadialMazeTurns - Compute distribution of successive turn angles in radial maze.
%
%	In order to detect stereotypic egocentric strategies that a rat could employ to
%  solve a radial maze task, this function determines the distribution (over trials)
%  of successive turn angles: if the rat does follow an egocentric strategy, then
%  successive turn angles are expected to reproduce a stereotypical sequence.
%
%  USAGE
%
%    distribution = RadialMazeTurns(data,mode)
%
%    data           list of [rat day trial arm configuration group] tuples
%    mode           optional complementary analysis; possibilities include
%                   'good' (include only trials with <=2 errors), 'preferred'
%                   (include only trials starting with the preferred arm) and
%                   'non-preferred' (include only trials starting with any
%                   arm except the preferred arm) (default = 'preferred')
%
%  NOTES
%
%    The data is assumed to be ordered, i.e. sorted by rat, configuration, day, and trial.
%
%    The 'preferred' arm is the most visited arm on the first visit.

% Copyright (C) 2009-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global RAT DAY TRIAL CONF GROUP ARM;

if nargin < 2,
	mode = 'preferred';
else
	mode = lower(mode);
end

nTrials = 3;
nArms = 8;
turnAngles = (-3:4)'*45;
armAngles = mod(360-(2:9)'*45,360)-180; % arm 1 is at 12:00, arm 3 at 15:00, etc.

% Column names for the 'data' matrix
RAT = 1;
DAY = 2;
TRIAL = 3;
CONF = 4;
GROUP = 5;
ARM = 6;
arms = data(:,ARM);
trials = data(:,TRIAL);
rats = data(:,RAT);
configurations = data(:,CONF);
days = data(:,DAY);
groups = data(:,GROUP);

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help RadialMazeTurns">RadialMazeTurns</a>'' for details).');
end


% Turn angles are differences in successive arm angles (numbers)
turns = diff(data(:,ARM));
% Keep in [1,8]
turns = mod(4-turns,8);
turns(turns==0) = 8;
turns = [0;turns];

% Differences in succesive arm angles are meaningless when they correspond to different rats
% (or configurations, or days, or trials): we will have to remove those cases, so we store this
% information for later use
change_rat = diff(data(:,RAT)) ~= 0;
change_configuration = diff(data(:,CONF)) ~= 0;
change_day = diff(data(:,DAY)) ~= 0;
change_trial = diff(data(:,TRIAL)) ~= 0;
change = change_rat|change_configuration|change_day|change_trial;
change = [true;change];
valid = ~change;

% Trials longer than this limit will be discarded from the final analysis
maxVisits = 5;

% Group days 1-5, 6-10 and 11-15
nBlocks = 3;
blocks = Bin(days,nBlocks);

% For each trial, label visits from 1 to N
newTrial = [change;false];
nVisits = diff([1;find(newTrial)]);
visits = ones(size(change));
visits(newTrial) = -nVisits+1;
visits = cumsum(visits);

P = [];
nGroups = max(groups);
for group = 1:nGroups,

	ratsInThisGroup = unique(rats(groups==group));

	nRats = length(ratsInThisGroup);
	for r = 1:nRats,
		thisRat = rats==ratsInThisGroup(r);
		fig = figure('position',[1285 52 1272 890]);
		c1 = uicontainer('position',[0 0 1 0.85]);
		c2 = uicontainer('position',[0 0.85 1 0.08]);
		c3 = uicontainer('position',[0 0.93 1 0.07]);

		nConfigurations = max(configurations);
		for configuration = 1:nConfigurations,

			thisConfiguration = configurations==configuration;

			for block = 1:nBlocks+1,  % iteration nBlocks+1 will provide more details about last block

				thisBlock = blocks==min([block nBlocks]);

				% Find most visited arm on first visit
				nVisitsInEachArm = hist(arms(thisRat&thisConfiguration&visits==1),1:8);
				[unused,preferredArm] = max(nVisitsInEachArm);

				% Compute distributions of turns (across trials and days) for successive visits
				T = turns(valid&thisRat&thisConfiguration&thisBlock);
				N = visits(valid&thisRat&thisConfiguration&thisBlock)-1; % turn number
				if block == nBlocks+1,
					if strcmp(mode,'preferred'),
						% Discard trials starting with non-preferred arm
						a = arms([valid(2:end,1);logical(1)]&thisRat&thisConfiguration&thisBlock);
						i = length(N);
						while i > 1,
							n = N(i);
							if a(i-n+1) ~= preferredArm,
								T(i-n+1:i)=[];
								N(i-n+1:i)=[];
								a(i-n+1:i)=[];
							end
							i = i - n;
						end
					elseif strcmp(mode,'non-preferred'),
						% Discard trials starting with preferred arm
						a = arms([valid(2:end,1);logical(1)]&thisRat&thisConfiguration&thisBlock);
						i = length(N);
						while i > 1,
							n = N(i);
							if a(i-n+1) == preferredArm,
								T(i-n+1:i)=[];
								N(i-n+1:i)=[];
								a(i-n+1:i)=[];
							end
							i = i - n;
						end
					else
						% Discard long trials
						i = length(N);
						while i > 1,
							n = N(i);
							if n > maxVisits-1,
								T(i-n+1:i)=[];
								N(i-n+1:i)=[];
							end
							i = i - n;
						end
					end
					if length(T) < 2, continue; end
				end
				turnDistribution = Accumulate([T N],1,[8 max(N)]);
				nTurns = size(turnDistribution,2);

				% Plot arm histograms
				A = arms(thisRat&thisConfiguration&thisBlock);
				V = visits(thisRat&thisConfiguration&thisBlock);
				if block == nBlocks+1,
					if strcmp(mode,'preferred'),
						% Discard trials starting with non-preferred arm
						i = length(V);
						a = arms([valid(1,2:end);logical(1)]&thisRat&thisConfiguration&thisBlock);
						while i > 1,
							n = V(i);
							if a(i-n+1) ~= preferredArm,
								V(i-n+1:i)=[];
								A(i-n+1:i)=[];
								a(i-n+1:i)=[];
							end
							i = i - n;
						end
					elseif strcmp(mode,'non-preferred'),
						% Discard trials starting with preferred arm
						a = arms([valid(2:end,1);logical(1)]&thisRat&thisConfiguration&thisBlock);
						i = length(N);
						while i > 1,
							n = N(i);
							if a(i-n+1) == preferredArm,
								V(i-n+1:i)=[];
								A(i-n+1:i)=[];
								a(i-n+1:i)=[];
							end
							i = i - n;
						end
					else
						% Discard long trials
						i = length(V);
						while i > 1,
							n = V(i);
							if n > maxVisits,
								V(i-n+1:i)=[];
								A(i-n+1:i)=[];
							end
							i = i - n;
						end
					end
				end
				armDistribution = Accumulate([A V],1,[8 max(V)]);
				nVisits = size(armDistribution,2);
				for i = 1:nVisits,
					subplot(1,nVisits*(nBlocks+1),i+(block-1)*nVisits,'parent',c3);hold on;
					%  Setup axes
					m = 1.5*max(armDistribution(:,i));
					tt = 0:0.01:2*pi;
					fill(m*cos(tt),m*sin(tt),'w');
					%  Plot
					h = rose(armAngles(A(V==i))*pi/180,(0:22.5:360)*pi/180);
					set(h,'color','k');
					xlim([-1 1]*m);
					ylim([-1 1]*m);
					axis square;axis off;
				end

				% Plot turn histograms
				for i = 1:nTurns,
					subplot(1,nTurns*(nBlocks+1),i+(block-1)*nTurns,'parent',c2);hold on;
					% Setup axes
					m = 1.5*max(turnDistribution(:,i));
					tt = 0:0.01:2*pi;
					fill(m*cos(tt),m*sin(tt),'w');
					% Plot
					h = rose(turnAngles(T(N==i))*pi/180,(0:22.5:360)*pi/180);
%  					if sum(N==i) >= 1,
%  						H = circ_raotest(turnAngles(T(N==i))*pi/180);
%  						if H == 1, set(h,'color','r'); end
%  					end
					p = circ_rtest(turnAngles(T(N==i))*pi/180);
					if p < 0.05, set(h,'color','r'); end
					if block == nBlocks+1,
						P = [P;p i group];
					end
					xlim([-1 1]*m);
					ylim([-1 1]*m);
					axis square;axis off;
				end

%  %  				subplot(2,4,block,'parent',c1);
%  %  				colormap(jet(6));
%  %  				td = [turnDistribution;turnDistribution(end,:)];td = [td td(:,end)];
%  %  				pcolor(0:size(turnDistribution,2),(-3:5)*45,td);
%  %  				caxis([0 6]);
%  %  				set(gca,'ytick',(-3:4)*45);
%  %  %  				set(gca,'xtick',[],'ytick',[]);
%  %  				shading('flat');

				% Plot distribution of turns
				subplot(2,4,block,'parent',c1);
				colormap(jet(6));
				image(1:size(turnDistribution,2),(-3:4)*45,turnDistribution+1);
				caxis([0 6]);set(gca,'ydir','normal','ytick',(-3:4)*45);

				% Plot successive turns (across trials and days)
				subplot(2,4,block,'parent',c1);
				hold on;
				start = [find(N==1);length(N)+1];
				for i = 1:length(start)-1,
					k = start(i):start(i+1)-1;
					plot(N(k),turnAngles(T(k))+rand(length(k),1)*15,'w','LineWidth',2);
				end

				% Plot turn(i+1) = f(turn(i))
				subplot(2,4,4+block,'parent',c1);
				hold on;
				start = [find(N==1);length(N)+1];
				for i = 1:length(start)-1,
					x = turnAngles(T(start(i):start(i+1)-2));
					y = turnAngles(T(start(i)+1:start(i+1)-1));
					plot(x+rand(length(x),1)*15,y+rand(length(x),1)*15,'k.-','LineWidth',1,'MarkerSize',12);
				end

				% Fix axes
				subplot(2,4,4+block,'parent',c1);xlim([-155 200]);ylim([-155 200]);
				set(gca,'xtick',turnAngles,'ytick',turnAngles);
			end
		end
		xlabel(['Group ' int2str(group) ' - Rat ' int2str(ratsInThisGroup(r)) ' [' mode ']']);
	end
end

figure;hold on;
colors = 'rbkgmy';
for group = 1:nGroups,
	p = P(P(:,3)==group,:);
	plot((p(:,3)-1)*10+p(:,2),p(:,1),colors(group),'marker','.','markersize',16,'linestyle','none');
end
xlim([0 nGroups*10]);
PlotHVLines(0.05,'h',':');