function Browse(parameter1,parameter2)

%Browse - Interactively browse the data in one or all subplots in a figure.
%
%  USAGE
%
%    Browse(figure,'on')   % enable browsing for this figure
%    Browse(figure,'off')  % disable browsing for this figure
%    Browse('on')          % enable browsing for all new figures
%    Browse('off')         % disable browsing for new figures
%    Browse('all')         % enable browsing for all existing figures
%    Browse('none')        % disable browsing for all existing figures
%
%    Browse(figure)        % same as Browse(figure,'on')
%    Browse                % same as Browse(gcf,'on')
%
%  NOTE
%
%    It is generally a good idea to append "Browse('on');" in the <a href="matlab:edit startup.m">startup</a> file,
%    so that it is automatically enabled for all figures as soon as Matlab is
%    started.
%
%  INTERACTIVE COMMANDS
%
%    To browse the data, press a key (hold SHIFT to browse all subplots):
%
%     left arrow     move one step to the left
%     right arrow    move one step to the right
%     page up        move one page to the left
%     page down      move one page to the right
%     home           move to the beginning
%     end            move to the end
%     i              zoom in (halve x axis width)
%     o              zoom out (double x axis width)
%
%    To graphically define new axis limits, click and drag the region of interest,
%    then press a key to update the axes (hold SHIFT to update all subplots):
%
%     x              change X axis
%     y              change Y axis
%     b              change both axes
%
%    To reset the axis limits, press a key (hold SHIFT to reset all subplots):
%
%     CTRL+x         reset X axis
%     CTRL+y         reset Y axis
%     CTRL+b         reset both axes
%
%    To use a special interactive commmand, click and drag the region of interest,
%    then press a key:
%
%       m            show mean and std of data in selected area (image only)
%       s            show X and Y coordinates of selected area

% Copyright (C) 2004-2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Name of this function
name = mfilename;
% Current default figure creation function
default = get(0,'DefaultFigureCreateFcn');
if isa(default,'function_handle'), default = func2str(default); end

% Special cases: Browse('on') and Browse('off')
if nargin == 1 && isstring_FMAT(lower(parameter1),'on','off','all','none'),
	% Current default figure creation function
	switch lower(parameter1),
		case 'on',
			if isempty(regexp(default,name)),
				set(0,'DefaultFigureCreateFcn',[name ';' default]);
			end
			help
		case 'off',
			if ~isempty(regexp(default,name)),
				set(0,'DefaultFigureCreateFcn',regexprep(default,[name ';'],''));
			end
		case 'all',
			children = get(0,'children');
			for i = 1:length(children),
				set(children(i),'WindowButtonDownFcn',@DefineSelection);
				set(children(i),'KeyPressFcn',@KeyPress);
			end
		case 'none',
			children = get(0,'children');
			for i = 1:length(children),
				set(children(i),'WindowButtonDownFcn',[]);
				set(children(i),'KeyPressFcn',[]);
			end
	end
	return
end

% Special case: call upon figure creation
if ~isempty(gcbo),
	% Current figure creation function
	current = get(gcbo,'CreateFcn');
	if isa(current,'function_handle'), current = func2str(current); end
	set(gcbo,'CreateFcn',regexprep(current,[name ';'],''));
end

% General case: parameter1 is a figure handle
if nargin == 0,
	fig = gcf;
else
	fig = parameter1;
end

if nargin == 2 && strcmp(lower(parameter2),'off')
	set(fig,'WindowButtonDownFcn',[]);
	set(fig,'KeyPressFcn',[]);
else
	set(fig,'WindowButtonDownFcn',@DefineSelection);
	set(fig,'KeyPressFcn',@KeyPress);
end

%-----------------------------------------------------------------------------------------------------------------
%   Help message
%-----------------------------------------------------------------------------------------------------------------

function help

disp(' ------------------------------------------------------------------------------');
disp('  Browsing mode automatically enabled for new figures.');
disp('  ');
disp('  To browse current axes, press a key (hold SHIFT to browse all subplots):');
disp('  ');
disp('     left arrow     move one step to the left');
disp('     right arrow    move one step to the right');
disp('     page up        move one page to the left');
disp('     page down      move one page to the right');
disp('     home           move to the beginning');
disp('     end            move to the end');
disp('     i              zoom in (halve x axis width)');
disp('     o              zoom out (double x axis width)');
disp('  ');
disp('  To graphically define new axis limits, click and drag the region of interest,');
disp('  then press a key to update the axes (hold SHIFT to update all subplots):');
disp('  ');
disp('     x              change X axis');
disp('     y              change Y axis');
disp('     b              change both axes');
disp('  ');
disp('  To reset the axis limits, press a key (hold SHIFT to reset all subplots):');
disp('  ');
disp('     CTRL+x         reset X axis');
disp('     CTRL+y         reset Y axis');
disp('     CTRL+b         reset both axes');
disp('  ');
disp('  To copy/paste axis limits, press a key:');
disp('  ');
disp('     c              copy current axis limits');
disp('     v              apply (paste) previously copied limits to current axis');
disp('  ');
disp('  To use a special interactive commmand, click and drag the region of interest,');
disp('  then press a key:');
disp('  ');
disp('     m               show mean and std over selected area (image only)');
disp('     s               show coordinates of selected area');
disp(' ------------------------------------------------------------------------------');

%-----------------------------------------------------------------------------------------------------------------
%   Mouse click and drag
%-----------------------------------------------------------------------------------------------------------------

function DefineSelection(fig,event)

if isempty(get(fig,'CurrentAxes')), return; end

% Let user select rectangle and store info until a key is pressed
point1 = get(gca,'CurrentPoint');
rbbox;
point2 = get(gca,'CurrentPoint');
x1 = point1(1,1);
x2 = point2(1,1);
xLim = [min([x1 x2]) max([x1 x2])];
y1 = point1(1,2);
y2 = point2(1,2);
yLim = [min([y1 y2]) max([y1 y2])];
set(fig,'UserData',[xLim;yLim]);

%-----------------------------------------------------------------------------------------------------------------
%   Key press
%-----------------------------------------------------------------------------------------------------------------

function KeyPress(fig,event)

if isempty(get(fig,'CurrentAxes')), return; end

if event.Key == 's',
	DisplaySelection(fig,event);
elseif event.Key == 'm',
	PlotSelectionMean(fig,event);
elseif event.Key == 'c',
	CopyAxisLimits(fig,event);
elseif event.Key == 'v',
	PasteAxisLimits(fig,event);
else
	ChangeLims(fig,event);
end

%-----------------------------------------------------------------------------------------------------------------
%   Display selection X and Y coordinates
%-----------------------------------------------------------------------------------------------------------------

function DisplaySelection(fig,event)

if isempty(get(fig,'CurrentAxes')), return; end

% Retrieve limits, based on previous mouse selection
lims = get(fig,'UserData');
if isempty(lims), return; end
xLim = lims(1,:);
yLim = lims(2,:);

disp(['X [' num2str(xLim(1)) '   ' num2str(xLim(2)) ']']);
disp(['Y [' num2str(yLim(1)) '   ' num2str(yLim(2)) ']']);

%-----------------------------------------------------------------------------------------------------------------
%   Plot mean and std of data enclosed in selection
%-----------------------------------------------------------------------------------------------------------------

function PlotSelectionMean(fig,event)

if isempty(get(fig,'CurrentAxes')), return; end

% Find image in gca (if any)
children = get(gca,'children');
img = [];
for i = 1:length(children),
	if strcmp(get(children(i),'type'),'image'),
		img = children(i);
	end
end
if isempty(img), return; end

% Retrieve limits, based on previous mouse selection
lims = get(fig,'UserData');
if isempty(lims), return; end
xLim = lims(1,:);
yLim = lims(2,:);
% Compute and plot mean
data = get(img,'cdata');
x = get(img,'xdata');
x = linspace(x(1),x(end),size(data,2));
y = get(img,'ydata');
y = linspace(y(1),y(end),size(data,1));
xOK = x >= xLim(1) & x <= xLim(2);
yOK = y >= yLim(1) & y <= yLim(2);
m = mean(data(yOK,xOK),2);
s = std(data(yOK,xOK),0,2);
figure;
PlotMean(y(yOK),m,m-s,m+s,':');

%-----------------------------------------------------------------------------------------------------------------
%   Copy current axis limits
%-----------------------------------------------------------------------------------------------------------------

function CopyAxisLimits(fig,event)

% (using weird names makes it unlikely that these globals could be used somewhere else)
global xAxisLimits_Browse49 yAxisLimits_Browse49;

xAxisLimits_Browse49 = xlim;
yAxisLimits_Browse49 = ylim;
disp(['Copy X axis [' num2str(xAxisLimits_Browse49(1)) '   ' num2str(xAxisLimits_Browse49(2)) ']']);
disp(['Copy Y axis [' num2str(yAxisLimits_Browse49(1)) '   ' num2str(yAxisLimits_Browse49(2)) ']']);

%-----------------------------------------------------------------------------------------------------------------
%   Apply (paste) previously copied limits to current axis
%-----------------------------------------------------------------------------------------------------------------

function PasteAxisLimits(fig,event)

% (using weird names makes it unlikely that these globals could be used somewhere else)
global xAxisLimits_Browse49 yAxisLimits_Browse49;

if ~isempty(xAxisLimits_Browse49) && ~isempty(yAxisLimits_Browse49),
	disp(['Paste X axis [' num2str(xAxisLimits_Browse49(1)) '   ' num2str(xAxisLimits_Browse49(2)) ']']);
	disp(['Paste Y axis [' num2str(yAxisLimits_Browse49(1)) '   ' num2str(yAxisLimits_Browse49(2)) ']']);
	xlim(xAxisLimits_Browse49);
	ylim(yAxisLimits_Browse49);
end

%-----------------------------------------------------------------------------------------------------------------
%   Change axes limits
%-----------------------------------------------------------------------------------------------------------------

function ChangeLims(fig,event)

if isempty(get(fig,'CurrentAxes')), return; end

% Define a few constants
NONE = 0;
RIGHT_PAGE = 1;
LEFT_PAGE = 2;
RIGHT_STEP = 3;
LEFT_STEP = 4;
HOME = 5;
END = 6;

% Default behavior (if no modifiers are pressed): modify only current axes, based on previous mouse selection
children = gca;
auto = false;
for i = 1:length(event.Modifier),
	if strcmp(event.Modifier{i},'shift'),
		% Shift pressed => change all subplots
		children = get(fig,'children');
	elseif strcmp(event.Modifier{i},'control'),
		% Control pressed => reset selected axes
		auto = true;
	end
end

% Which action?
changeX = false;
changeY = false;
move = NONE;
zoomX = 1;
switch(event.Key),
	case 'x',
		changeX = true;
	case 'y',
		changeY = true;
	case 'b',
		changeX = true;
		changeY = true;
	case 'pagedown',
		move = RIGHT_PAGE;
	case 'pageup',
		move = LEFT_PAGE;
	case 'rightarrow',
		move = RIGHT_STEP;
	case 'leftarrow',
		move = LEFT_STEP;
	case 'home',
		move = HOME;
	case 'end',
		move = END;
	case 'i',
		zoomX = 0.5;
	case 'o',
		zoomX = 2;
	otherwise,
		return
end

% Update subplot(s)
for i = 1:length(children),
	type = get(children(i),'type');
	if strcmp(type,'axes'),
		if move ~= NONE,
			% Move x lims
			x = xlim(children(i));
			dx = x(2)-x(1);
			switch(move),
				case RIGHT_PAGE,
					xLim = x+dx;
				case LEFT_PAGE,
					xLim = x-dx;
				case RIGHT_STEP,
					xLim = x+dx/10;
				case LEFT_STEP,
					xLim = x-dx/10;
				case HOME,
					% Determine min x value for all data
					c = get(children(i),'children');
					minX = Inf;
					for j = 1:length(c),
						m = get(c(j),'xdata');
						minX = min([minX m(1)]);
					end
					xLim = minX + [0 dx];
				case END,
					% Determine max x value for all data
					c = get(children(i),'children');
					maxX = -Inf;
					for j = 1:length(c),
						m = get(c(j),'xdata');
						maxX = max([maxX m(end)]);
					end
					xLim = maxX - [dx 0];
			end
			set(children(i),'xlim',xLim);
			% Reset mouse selection area
			set(fig,'UserData',[]);
		elseif zoomX ~= 1,
			% Zoom x axis
			xLim = get(children(i),'xlim');
			set(children(i),'xlim',[xLim(1) xLim(1)+zoomX*(xLim(2)-xLim(1))]);
		else
			% Retrieve limits, based on previous mouse selection
			lims = get(fig,'UserData');
			if ~auto && isempty(lims), return; end
			if ~isempty(lims),
				xLim = lims(1,:);
				yLim = lims(2,:);
			end
			% Change x and/or y lims
			if changeX,
				if auto,
					set(children(i),'xlimmode','auto');
				else
					if xLim(1) ~= xLim(2), set(children(i),'xlim',xLim); end
				end
			end
			if changeY,
				if auto,
					set(children(i),'ylimmode','auto');
				else
					if yLim(1) ~= yLim(2), set(children(i),'ylim',yLim); end
				end
			end
		end
	end
end