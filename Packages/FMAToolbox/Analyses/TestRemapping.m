function [h,p,cr,ca] = TestRemapping(control,repeat,test,varargin)

%TestRemapping - Test if firing fields remap (or shift) between two conditions.
%
%  To test for random remapping between a control and a test condition, we
%  first estimate for each place cell the absolute field shift between the two
%  conditions, as the mode of the spatial cross-correlogram between the
%  respective firing fields. We then divide this value by the average field
%  size across conditions, yielding a relative shift. This measures by how much
%  the field moves between the two conditions, in proportion to the field size.
%  To account for baseline variability, we also estimate the relative field
%  shift between two repetitions of the control condition, and subtract it from
%  the relative shift in the test condition. This (unsigned) corrected shift is
%  averaged over all cells.
%
%  To test if this is significantly different from random remapping, we generate
%  a distribution of corrected relative shifts under the null hypothesis of
%  random remapping. This is done by randomly permuting the test firing fields
%  between the cells (bootstrap).
%
%  The null hypothesis is rejected if the probability of the observed shift is
%  lower than the critical value.
%
%  In addition, we compute the (1-alpha) confidence interval of the (signed)
%  corrected shifts. These can be used to test for systematic (expected) field
%  shifts.
%
%  Current implementation only tests 1D environments.
%
%  USAGE
%
%    [h,p,cr,ca] = TestRemapping(control,repeat,test,<options>)
%
%    control        firing fields in control condition (MxN: M fields, N bins)
%    repeat         firing fields in repeated control condition
%    test           firing fields in test condition
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'type'        'linear' for linear tracks (default), 'circular' otherwise
%     'alpha'       significance level (default = 0.05)
%     'iterations'  number of iterations for bootstrap (default = 150)
%     'show'        plot results (default = 'off')
%    =========================================================================
%
%  OUTPUT
%
%    h              0 if H0 (random remapping) cannot be rejected, 1 otherwise
%    p              p-value of bootstrap test
%    cr             confidence interval for relative shifts
%    ca             confidence interval for absolute shifts

% Copyright (C) 2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
type = 'linear';
show = 'off';
nBootstrap = 150;
alpha = 0.05;

% Check number of parameters
if nargin < 3 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help TestRemapping">TestRemapping</a>'' for details).');
end

% Check parameter sizes
if ~isdmatrix(control) || ~isdmatrix(repeat) || ~isdmatrix(test),
	error('All firing fields should be MxN matrices (type ''help <a href="matlab:help TestRemapping">TestRemapping</a>'' for details).');
end
if ~(all(size(control)==size(repeat)) && all(size(control)==size(test))),
	error('All firing fields should have the same sizes (type ''help <a href="matlab:help TestRemapping">TestRemapping</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help TestRemapping">TestRemapping</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'type',
			type = varargin{i+1};
			if ~isstring_FMAT(type,'linear','circular'),
				error('Incorrect value for property ''type'' (type ''help <a href="matlab:help TestRemapping">TestRemapping</a>'' for details).');
			end
		case 'alpha',
			alpha = varargin{i+1};
			if ~isdscalar(alpha,'>=0','<=1'),
				error('Incorrect value for property ''alpha'' (type ''help <a href="matlab:help TestRemapping">TestRemapping</a>'' for details).');
			end
		case 'iterations',
			nBootstrap = varargin{i+1};
			if ~isiscalar(nBootstrap,'>0'),
				error('Incorrect value for property ''iterations'' (type ''help <a href="matlab:help TestRemapping">TestRemapping</a>'' for details).');
			end
		case 'show',
			show = varargin{i+1};
			if ~isstring_FMAT(show,'on','off'),
				error('Incorrect value for property ''show'' (type ''help <a href="matlab:help TestRemapping">TestRemapping</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help TestRemapping">TestRemapping</a>'' for details).']);
	end
end
[m,n] = size(control);

% Compute relative field shift (= relative to size)
[relativeShiftControl,absoluteShiftControl,xcControl] = FieldShift(control,repeat,'type',type);
[relativeShiftTest,absoluteShiftTest,xcTest] = FieldShift(control,test,'type',type);

% Mean relative shift, taking into account baseline variability
meanUnsignedRelativeShift = mean(abs(relativeShiftTest));

% Cumulative density function (CDF) of null distribution, estimated by bootstrap
global RemappingTest_control RemappingTest_shiftControl RemappingTest_type;
RemappingTest_control = control;
RemappingTest_shiftControl = relativeShiftControl;
RemappingTest_type = type;
bRelativeShift = bootstrp(nBootstrap,@RemappingShift,test); % see RemappingShift below
[bPDF,x] = hist(bRelativeShift,100);
bPDF = bPDF / nBootstrap;
bCDF = cumsum(bPDF);

% Test for random remapping
if meanUnsignedRelativeShift < min(x),
	p = 0;
elseif meanUnsignedRelativeShift > max(x),
	p = 1;
else
	p = interp1(x,bCDF,meanUnsignedRelativeShift);
end
h = double(p<alpha);

% Confidence intervals for relative and absolute shifts (via bootstrap)
s = bootstrp(nBootstrap,@mean,absoluteShiftTest);
ca = [prctile(s,100*(alpha/2)) prctile(s,100*(1-(alpha/2)))];
s = bootstrp(nBootstrap,@mean,relativeShiftTest);
cr = [prctile(s,100*(alpha/2)) prctile(s,100*(1-(alpha/2)))];

% Find peak positions (optimized code, uneasy to read)
[y,x] = find(control==repmat(max(control,[],2),1,n));
peakControl = Accumulate(y,x)./Accumulate(y,1)/n;
[y,x] = find(test==repmat(max(test,[],2),1,n));
peakTest = Accumulate(y,x)./Accumulate(y,1)/n;

if strcmp(show,'on'),

	figure;
	% Plot control cross-correlograms
	subplot(1,2,1);hold on;
	PlotColorCurves(xcControl,[-1 1],0,'w--',absoluteShiftControl,'.k','gamma',2);
	title(['Control (N=' int2str(size(xcControl,1)) ')']);
	% Plot test cross-correlograms
	subplot(1,2,2);hold on;
	PlotColorCurves(xcTest,[-1 1],0,'w--',absoluteShiftTest,'.k','gamma',2);
	title(['Test (N=' int2str(size(xcTest,1)) ')']);

   figure;
	% Plot control cross-correlograms in ascending order
	[absoluteShiftControl,order] = sort(absoluteShiftControl);
	xcControl = xcControl(order,:);
	subplot(1,2,1);hold on;
	PlotColorCurves(xcControl,[-1 1],0,'w--',absoluteShiftControl,'.k','gamma',2);
	% Plot histogram
	SideAxes('top',0.1,'gap',0.01);
	[hi,xh] = hist(absoluteShiftControl,50);
	bar(xh,hi);
	PlotHVLines(mean(absoluteShiftControl),'v','r','linewidth',2);
	xlim([-1 1]);
	set(gca,'xtick',[]);
	title(['Control (N=' int2str(size(xcControl,1)) ')']);
	% Plot test cross-correlograms in ascending order
	[absoluteShiftTest,order] = sort(absoluteShiftTest);
	xcTest = xcTest(order,:);
	subplot(1,2,2);hold on;
	PlotColorCurves(xcTest,[-1 1],0,'w--',absoluteShiftTest,'.k','gamma',2);
	% Plot histogram
	SideAxes('top',0.1,'gap',0.01);
	[hi,xh] = hist(absoluteShiftTest,50);
	bar(xh,hi);
	PlotHVLines(mean(absoluteShiftTest),'v','r','linewidth',2);
	xlim([-1 1]);
	set(gca,'xtick',[]);
	title(['Test (N=' int2str(size(xcTest,1)) ')']);

	figure;hold on;
	% Plot relative shift histogram
	[hi,xh] = hist(relativeShiftTest,50);
	bar(xh,hi/nBootstrap);
	PlotHVLines(cr,'v','k--');
	PlotHVLines(mean(s),'v','r','linewidth',2);
	xlabel('Relative Shift');
	ylabel('Probability');
	title(['N = ' num2str(nBootstrap)]);
	xLim = xlim;
	yLim = ylim;
	[~,pval] = ttest(relativeShiftTest);
	prop1 = sum(abs(relativeShiftTest)<=1)/length(relativeShiftTest);
	prop2 = sum(abs(relativeShiftTest)<=0.5)/length(relativeShiftTest);
	m = mean(relativeShiftTest);
	s = sem(relativeShiftTest);
	l = length(relativeShiftTest);
	text(xLim(1)+0.1*abs(xLim(1)),yLim(2)*0.9,{mean2str(m,s,l),[num2str((1-alpha)*100) '% Confidence Interval [' sprintf('%.3f',cr(1)) ',' num2str(cr(2)) ']'],['t-test p = ' sprintf('%.3f',pval)],[sprintf('%.3f',prop1*100) '% |x|<=1'],[sprintf('%.3f',prop2*100) '% |x|<=0.5']});

	figure;
	% Plot test vs control shifts (density plot)
	subplot(2,2,1);hold on;
	nBins = 500;
	range = max([relativeShiftControl;relativeShiftTest])*[-1 1];
	sc = Bin(relativeShiftControl,range,nBins);
	st = Bin(relativeShiftTest,range,nBins);
	z = Smooth(Accumulate([sc st],1,[nBins nBins]),10)'/n;
	PlotColorMap(z,'x',linspace(range(1),range(2),n),'y',linspace(range(1),range(2),n));
	xlabel('Control Shift (proportion of field size)');
	ylabel('Test Shift (proportion of field size)');
	clim([0 0.3*max(z(:))]);
	% Plot field modes in control vs test (density plot)
	subplot(2,2,2);hold on;
	PlotRepeat(peakControl,peakTest,'xxyy','k.');
	xlabel('Control Peak Position');
	ylabel('Test Peak Position');
	[beta,R2,pSlope] = CircularRegression(peakControl*2*pi,peakTest*2*pi,'slope',1);
	title(['slope = ' num2str(beta(1)) ' (p = ' num2str(pSlope) ')']);
	PlotSlope(1,1,beta(1),2,'r','LineWidth',2);
	% Plot shift distribution under H0 (bootstrap)
	subplot(2,2,[3 4]);hold on;
	[hi,x] = hist(bRelativeShift,30);
	b = bar(x,hi/nBootstrap);
	hold on;
	PlotHVLines(meanUnsignedRelativeShift,'v','r');
	Y = ylim; Y = Y(2);
	text(meanUnsignedRelativeShift*1.25,Y*0.85,['p = ' sprintf('%.3f',p)]);
	xlim([0 1]);
	title(['N = ' num2str(nBootstrap)]);
	xlabel('Mean Unsigned Shift (proportion of field size)');
	ylabel('Probability');
	
	figure;
	% Plot shift as a function of peak position
	subplot(2,2,1);hold on;
	pc = Bin(peakControl,[0 1],10);
	[m,stdev] = Accumulate(pc,relativeShiftTest);
	s = stdev/sqrt(length(m));
	%errorbar(0.05+(0:0.1:0.9),m,ci(:,1)-m,ci(:,2)-m,'k','LineStyle','none');
	errorbar(0.05+(0:0.1:0.9),m,s,s,'k','LineStyle','none');
	bar(0.05+(0:0.1:0.9),m);
	plot(peakControl,relativeShiftTest,'.r');
	xlabel('Peak Position');
	ylabel('Shift (proportion of field size)');
	pa = anova1(relativeShiftTest,pc,'off');
	yLim = ylim;
	text(mean(xlim),yLim(2)*0.9,['p = ' num2str(pa)]);
	subplot(2,2,3);hold on;
	[m,stdev] = Accumulate(pc,abs(relativeShiftTest));
	s = stdev/sqrt(length(m));
	errorbar(0.05+(0:0.1:0.9),m,s,s,'k','LineStyle','none');
	bar(0.05+(0:0.1:0.9),m);
	plot(peakControl,relativeShiftTest,'.r');
	xlabel('Peak Position');
	ylabel('Unsigned Shift (proportion of field size)');
	pa = anova1(abs(relativeShiftTest),pc,'off');
	yLim = ylim;
	text(mean(xlim),yLim(2)*0.9,['p = ' num2str(pa)]);
	subplot(2,2,2);hold on;
	[m,stdev] = Accumulate(pc,absoluteShiftTest);
	s = stdev/sqrt(length(m));
	errorbar(0.05+(0:0.1:0.9),m,s,s,'k','LineStyle','none');
	bar(0.05+(0:0.1:0.9),m);
	plot(peakControl,absoluteShiftTest,'.r');
	xlabel('Peak Position');
	ylabel('Shift (proportion of field size)');
	pa = anova1(absoluteShiftTest,pc,'off');
	yLim = ylim;
	text(mean(xlim),yLim(2)*0.9,['p = ' num2str(pa)]);
	subplot(2,2,4);hold on;
	[m,stdev] = Accumulate(pc,abs(absoluteShiftTest));
	s = stdev/sqrt(length(m));
	errorbar(0.05+(0:0.1:0.9),m,s,s,'k','LineStyle','none');
	bar(0.05+(0:0.1:0.9),m);
	plot(peakControl,absoluteShiftTest,'.r');
	xlabel('Peak Position');
	ylabel('Unsigned Shift (proportion of field size)');
	pa = anova1(abs(absoluteShiftTest),pc,'off');
	yLim = ylim;
	text(mean(xlim),yLim(2)*0.9,['p = ' num2str(pa)]);

end

% ------------------------------- Helper function -------------------------------

function meanUnsignedRelativeShift = RemappingShift(test)

global RemappingTest_control RemappingTest_shiftControl RemappingTest_type;
control = RemappingTest_control;
relativeShiftControl = RemappingTest_shiftControl;
type = RemappingTest_type;

% Compute field shift as the distance relative to size
relativeShiftTest = FieldShift(control,test,'type',type);

% Shift, taking into account 'intrinsic' variability
meanUnsignedRelativeShift = mean(abs(relativeShiftTest));
