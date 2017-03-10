function status = Hide(parameter1,parameter2)

%Hide - Hide (or show) existing or future figures, e.g. to speed up batch processing.
%
%  USAGE
%
%    Hide(figures,'on')   % hide these figures
%    Hide(figures,'off')  % un-hide this figure
%    Hide('all')          % hide all existing figures
%    Hide('none')         % un-hide all existing figures
%    Hide('on')           % automatically hide all new figures
%    Hide('off')          % do not automatically hide new figures
%
%    Hide(figures)        % same as Hide(figures,'on')
%    Hide                 % same as Hide(gcf,'on')
%    Hide('status')       % returns the current default behavior ('on' or 'off')
%
%  SEE ALSO
%
%    See also sca, scf, StartBatch.
%

% Copyright (C) 2011-2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

status = [];

% Name of this function
name = mfilename;
% Current default figure creation function
default = get(0,'DefaultFigureCreateFcn');
if isa(default,'function_handle'), default = func2str(default); end

% Special cases: Hide('on'), Hide('off'), Hide('all'), Hide('none')
if nargin == 1 && isstring(lower(parameter1),'on','off','all','none','status'),
	switch lower(parameter1),
		case 'status',
			if isempty(regexp(default,name)),
				status = 'off';
			else
				status = 'on';
			end
		case 'on',
			if isempty(regexp(default,name)),
				set(0,'DefaultFigureCreateFcn',[name ';' default]);
			end
		case 'off',
			if ~isempty(regexp(default,name)),
				set(0,'DefaultFigureCreateFcn',regexprep(default,[name ';'],''));
			end
		case 'all',
			children = get(0,'children');
			for i = 1:length(children),
				if strcmp(get(children(i),'visible'),'on'),
					% Do not hide if already hidden, because this is *slow*
					set(children(i),'visible','off');
				end
			end
		case 'none',
			children = get(0,'children');
			for i = 1:length(children),
				if strcmp(get(children(i),'visible'),'off'),
					% Do not unhide if already visible, because this is *slow*
					set(children(i),'visible','on');
				end
			end
	end
	return
end

% Special case: call upon figure creation
if ~isempty(gcbo),
	% Make invisible
	set(gcbo,'visible','off');
	warning('FMAToolbox:Hide:FigureHidden','Figure hidden (type ''help <a href="matlab:help Hide">Hide</a>'' for details).');
	% Remove from figure creation function
	current = get(gcbo,'CreateFcn');
	if isa(current,'function_handle'), current = func2str(current); end
	set(gcbo,'CreateFcn',regexprep(current,[name ';'],''));
end

% General case: parameter1 is a list of figure handles
if nargin < 1,
	figs = gcf;
else
	figs = parameter1;
end
if nargin < 2,
	mode = 'on';
else
	mode = lower(parameter2);
end

switch mode,
	case 'off',
		for fig = figs,
			set(fig,'visible','on');
		end
	case 'on',
		for fig = figs,
			set(fig,'visible','off');
		end
end

