function stats = FieldPSP(stims,lfp,varargin)

%FieldPSP - Measure field EPSPs (waveforms, slope, amplitude) across time.
%
%  Slope is determined as the maximum slope during the time window around stimulation.
%  Amplitude is determined as the maximum amplitude during the same time window.
%
%  USAGE
%
%    stats = FieldPSP(stims,lfp,<options>)
%
%    stims          stimulation timestamps
%    lfp            local field potential samples
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'direction'   direction of the field EPSP, either 'up' (recording
%                   from str. pyr., default) or 'down' (recording from str.
%                   radiatum)
%     'durations'   durations before and after synchronizing events for each
%                   trial (in s) (default = [-0.01 0.02])
%     'slope'       where max slope should be sought (default = [0.001 0.005])
%     'amplitude'   where max amplitude should be sought
%                   (default = [0.005 0.010])
%     'mode'        compute slope and amplitude running averages using a
%                   fixed time window ('time'), or a window with a fixed
%                   number of stimulations ('count', default)
%     'window'      window length (in s or in counts) for running averages
%                   (default = 100)
%     'show'        either 'on' or 'off' (default)
%    =========================================================================
%
%  OUTPUT
%
%    stats.amplitude.t         times of max amplitude (relative to stimulation)
%    stats.amplitude.v         values of max amplitudes
%    stats.slope.t             times of max slopes (relative to stimulation)
%    stats.slope.v             values of max slopes
%    stats.psp.t               times of fPSP samples (relative to stimulation)
%    stats.psp.v               fPSP waveforms, one per stimulation
%

% Copyright (C) 2010-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
parent = [];
durations = [-0.01 0.02];
mode = 'count';
window = 100;
direction = 'up';
stats = [];
show = 'off';
slopeWindow = [0.001 0.005];
amplitudeWindow = [0.005 0.010];

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help FieldPSP">FieldPSP</a>'' for details).');
end

% Check parameter sizes
if ~isdvector(stims),
	error('Parameter ''stims'' is not a vector (type ''help <a href="matlab:help FieldPSP">FieldPSP</a>'' for details).');
end
if ~isdmatrix(lfp),
	error('Parameter ''lfp'' is not a matrix (type ''help <a href="matlab:help FieldPSP">FieldPSP</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FieldPSP">FieldPSP</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'direction',
			direction = lower(varargin{i+1});
			if ~isstring_FMAT(direction,'up','down'),
				error('Incorrect value for property ''direction'' (type ''help <a href="matlab:help FieldPSP">FieldPSP</a>'' for details).');
			end
		case 'durations',
			durations = varargin{i+1};
			if ~isdvector(durations,'#2','<'),
				error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help FieldPSP">FieldPSP</a>'' for details).');
			end
		case 'slope',
			slopeWindow = varargin{i+1};
			if ~isdvector(slopeWindow,'#2','<'),
				error('Incorrect value for property ''slope'' (type ''help <a href="matlab:help FieldPSP">FieldPSP</a>'' for details).');
			end
		case 'amplitude',
			amplitudeWindow = varargin{i+1};
			if ~isdvector(amplitudeWindow,'#2','<'),
				error('Incorrect value for property ''amplitude'' (type ''help <a href="matlab:help FieldPSP">FieldPSP</a>'' for details).');
			end
		case 'mode',
			mode = lower(varargin{i+1});
			if ~isstring_FMAT(mode,'time','count'),
				error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help FieldPSP">FieldPSP</a>'' for details).');
			end
		case 'parent',
			parent = varargin{i+1};
			if ~ishandle(parent),
				error('Incorrect value for property ''parent'' (type ''help <a href="matlab:help FieldPSP">FieldPSP</a>'' for details).');
			end
		case 'window',
			window = varargin{i+1};
			if ~isdscalar(window,'>0'),
				error('Incorrect value for property ''window'' (type ''help <a href="matlab:help FieldPSP">FieldPSP</a>'' for details).');
			end
		case 'show',
			show = lower(varargin{i+1});
			if ~isstring_FMAT(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help FieldPSP">FieldPSP</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FieldPSP">FieldPSP</a>'' for details).']);
	end
end

up = strcmp(direction,'up');
show = strcmp(show,'on');

% Some analyses depend on the method (fixed time vs fixed count), but others do not
% To simplify the code, we will use two variables:
%  - stims0 uses time, whatever the method
%  - stims  uses either time or count, depending on the method
stims0 = stims;
if strcmp(mode,'count'),
	stims = (1:length(stims))';
end

% Constants
nStims = length(stims);
height = max(lfp(:,2))-min(lfp(:,2));
nSlices = 12;
start = stims(1);
stop = stims(end);
window = (stop-start)/nSlices;
slices = (0:nSlices-1)*window+start;
nStims = length(stims);
frequency = 1/median(diff(lfp(:,1)));
smooth = frequency/1250;
nStimsPlottedPerSlice = 5;
nStimsPerSlice = nStims/nSlices;
nSkip = floor(nStims/(nSlices*(nStimsPlottedPerSlice-1))); % plot one stim every n
colors = Bright(nStimsPlottedPerSlice);

% Compute slope and amplitude for each stimulation

% Resynchronize all fPSPs (use only 15ms after stims)
[sync,index] = Sync(lfp,stims0,'durations',durations);

% Preallocate large matrix to store all fPSPs
n = sum(index==1);
stats.psp.v = nan(nStims,n);
stats.psp.t = linspace(durations(1),durations(2),n);
% Compute slope and amplitude for each stimulation
if show, figure; end
for i = 1:nStims,
	% Select the appropriate data, store in fPSP matrix, smooth and differentiate
	s = sync(index==i,:);
	stats.psp.v(i,:) = interp1(s(:,1),s(:,2),stats.psp.t);
	s(:,2) = Smooth(s(:,2),smooth);
	% Compute amplitude as the first local maximum (or minimum) in the appropriate time interval
	in = InIntervals(s(:,1),amplitudeWindow);
	if up,
		extrema = find(IsExtremum(s)&in);
	else
		extrema = find(IsExtremum(s,'mode','minima')&in);
	end
	if isempty(extrema),
		amplitude(i) = nan;
		tAmplitude(i) = nan;
		slope(i) = nan;
		tSlope(i) = nan;
		continue;
	end
	extremum = extrema(1);
	tAmplitude(i) = s(extremum,1);
	amplitude(i) = s(extremum,2);
	% Compute slope as the maximum slope in the appropriate time interval
	% (set all values outside this time interval to nan to make sure the max will be in the interval)
	ds = Diff(s);
	ds(~InIntervals(ds(:,1),slopeWindow),2) = nan;
	if up,
		[slope(i),j] = nanmax(ds(:,2));
	else
		[slope(i),j] = nanmin(ds(:,2));
	end
	tSlope(i) = ds(j,1);
	% Plot slope and amplitude
	slice = floor(i/nStimsPerSlice);
	stimInSlice = floor(i-1-slice*nStimsPerSlice);
	if show && rem(stimInSlice,nSkip) == 0 && stimInSlice > 0 && stimInSlice < nStimsPlottedPerSlice,
		k = ceil(nSlices*i/nStims);
		l = floor(stimInSlice/nSkip)+1;
		SquareSubplot(nSlices,k);
		xlabel(num2str(slices(k)));
		hold on;
		PlotXY(s,'color',colors(l,:),'linewidth',3);
		PlotSlope(tAmplitude(i),amplitude(i),0,0.005,'r');
		PlotSlope(tSlope(i),s(j,2),slope(i),0.005,'k');
	end
end

% Make sure all subplots use the same scale
if show,
	sub = get(gcf,'children');
	if ~isempty(sub),
		for i = 1:length(sub),
			lims(i,:) = ylim(sub(i));
		end
		m = min(lims(:,1));
		M = max(lims(:,2));
		for i = 1:length(sub),
			ylim(sub(i),[m M]);
		end
	end
end

stats.amplitude.t = tAmplitude;
stats.amplitude.v = amplitude;
stats.slope.t = tSlope;
stats.slope.v = slope;

% Additional figures: mean fPSPs, slopes and amplitudes (scatterplots + running averages), CVs, ISIs and short time CCGs

if show,
	% Mean fPSPs
	figure;
	for i = 1:nSlices,
		if strcmp(mode,'count'),
			start = floor((i-1)*nStimsPerSlice)+1;
			stop = floor(i*nStimsPerSlice);
			ok = start:stop;
		else
			ok = find(InIntervals(stims,[0 window]+start+(i-1)*window));
		end
		m = mean(stats.psp.v(ok,:));
		e = sem(stats.psp.v(ok,:));
		SquareSubplot(nSlices,i);
		xlabel(num2str(slices(i)));
		hold on;
		PlotMean(stats.psp.t,m,m-e/2,m+e/2,':');
	end

	figure;
	% 'Rectify' slope and amplitude
	if ~up,
		amplitude = -amplitude;
		slope = -slope;
	end
	if strcmp(mode,'time'),
		xLabel = 'Time (s)';
	else
		xLabel = 'Count';
	end
	% Slope (scatterplot)
	subplot(2,3,1);
	hold on;
	plot(stims,slope,'.','markersize',1);
	% Slope (running average)
	[ts,ms,es] = RunningAverage(stims,slope,'window',window,'overlap',0.5*window);
	PlotMean(ts,ms,es(:,1),es(:,2),':','k');
	PlotHVLines(slices,'v','color',[0.8 0.8 0.8]);
	xlim([stims(1,1) stims(end,1)]);
	ylabel('Slope (a.u.)');
	set(gca,'xtick',[]);
	% Slope time (scatterplot)
	subplot(2,3,2);
	hold on;
	plot(stims,tSlope*1000,'.','markersize',1);
	% Slope time (running average)
	[ts,ms,es] = RunningAverage(stims,tSlope*1000,'window',window,'overlap',0.5*window);
	PlotMean(ts,ms,es(:,1),es(:,2),':','k');
	PlotHVLines(slices,'v','color',[0.8 0.8 0.8]);
	xlim([stims(1,1) stims(end,1)]);
	ylabel('Slope time (ms)');
	set(gca,'xtick',[]);

	% Amplitude (scatter)
	subplot(2,3,4);
	hold on;
	plot(stims,amplitude,'.','markersize',1);
	% Amplitude (running average)
	[ta,ma,ea] = RunningAverage(stims,amplitude,'window',window,'overlap',0.5*window);
	PlotMean(ta,ma,ea(:,1),ea(:,2),':','k');
	PlotHVLines(slices,'v','color',[0.8 0.8 0.8]);
	xlim([stims(1,1) stims(end,1)]);
	xlabel(xLabel);
	ylabel('Amplitude (a.u.)');
	% Amplitude time (scatter)
	subplot(2,3,5);
	hold on;
	plot(stims,tAmplitude*1000,'.','markersize',1);
	% Amplitude time (running average)
	[ta,ma,ea] = RunningAverage(stims,tAmplitude*1000,'window',window,'overlap',0.5*window);
	xlim([stims(1,1) stims(end,1)]);
	PlotMean(ta,ma,ea(:,1),ea(:,2),':','k');
	PlotHVLines(slices,'v','color',[0.8 0.8 0.8]);
	ylabel('Amplitude time (ms)');
	xlabel(xLabel);

	% ISIs
	subplot(2,3,3);
	ds = diff(stims0(:,1));
	x = 0:0.1:10;
	h = hist(ds,x);
	h = h/sum(h);
	bar(x,h);
	xlabel('ISI');
	xlim([0 10]);
	ylabel('Count');
	% CV and CV2
	binSize = 0.1;
	smooth = 10;
	for i = 1:20,
		cv(i) = CV(stims0,'order',i,'binSize',binSize,'smooth',smooth,'method','fixed');
	end
	cv2 = CV(stims0,'measure','cv2');
	subplot(2,3,6);
	plot(cv,'+-');
	hold on;
	plot(cv2,'r+');
	xlabel('nth ISI');
	ylabel('CV & CV2');

	% Stims autocorrelograms
	figure;
	[ccg,x,y] = ShortTimeCCG(stims0,'min',1,'mode','norm','smooth',[2 0]);
	PlotShortTimeCCG(ccg,'x',x,'y',y);
	clim([0 0.01]);
end