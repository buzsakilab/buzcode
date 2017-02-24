function AdjustAxes(f,varargin)

%AdjustAxes - Adjust axes limits for all subplots
%
%  USAGE
%
%    AdjustAxes(f,<options>)
%
%    f              optional figure handle (default = gcf)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'x'           [min max], 'uniform' (same for all axes) or 'auto'
%     'y'           [min max], 'uniform' (same for all axes) or 'auto'
%     'xy'          [min max], 'uniform' (same for all axes) or 'auto'
%    =========================================================================
%
%  DEFAULTS
%
%    When no option is provided, all axes are automatically adjusted (same as
%    xlim('xy','uniform')). When 'x' is provided alone, y axes are left unchanged
%    and vice versa.
%

% Copyright (C) 2012 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
xlims = [];
ylims = [];

% Optional parameter
if nargin < 1,
	f = gcf;
elseif ischar(f),
	varargin = {f,varargin{:}};
	f = gcf;
elseif ~ishandle(f),
	error('Incorrect value for ''f'' (type ''help <a href="matlab:help AdjustAxes">AdjustAxes</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
  if ~ischar(varargin{i}),
    error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help AdjustAxes">AdjustAxes</a>'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'x',
      xlims = lower(varargin{i+1});
      if ~isdvector(xlims,'<') && ~isstring_FMAT(xlims,'auto','uniform'),
        error('Incorrect value for property ''xlims'' (type ''help <a href="matlab:help AdjustAxes">AdjustAxes</a>'' for details).');
      end
    case 'y',
      ylims = lower(varargin{i+1});
      if ~isdvector(ylims,'<') && ~isstring_FMAT(ylims,'auto','uniform'),
        error('Incorrect value for property ''ylims'' (type ''help <a href="matlab:help AdjustAxes">AdjustAxes</a>'' for details).');
      end
    case 'xy',
      xlims = lower(varargin{i+1});
      ylims = lower(varargin{i+1});
      if ~isdvector(xlims,'<') && ~isstring_FMAT(xlims,'auto','uniform'),
        error('Incorrect value for property ''AdjustAxess'' (type ''help <a href="matlab:help AdjustAxes">AdjustAxes</a>'' for details).');
      end
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help AdjustAxes">AdjustAxes</a>'' for details).']);
  end
end

if isempty(xlims) && isempty(ylims),
	xlims = 'uniform';
	ylims = 'uniform';
end

% Compute min and max limits for all subplots
sub = get(f,'children');
if strcmp(xlims,'uniform'),
	lims = [];
	for i = 1:length(sub),
		if strcmp(get(sub,'type'),'axes'),
			lims(end+1,:) = xlim(sub(i));
		end
	end
	if isempty(lims), return; end
	xlims = [min(lims(:,1)) max(lims(:,2))];
end
if strcmp(ylims,'uniform'),
	lims = [];
	for i = 1:length(sub),
		if strcmp(get(sub,'type'),'axes'),
			lims(end+1,:) = ylim(sub(i));
		end
	end
	if isempty(lims), return; end
	ylims = [min(lims(:,1)) max(lims(:,2))];
end

% Apply
for i = 1:length(sub),
	if strcmp(get(sub,'type'),'axes'),
		if strcmp(xlims,'auto'),
			xlim(sub(i),'auto');
		elseif ~isempty(xlims),
			xlim(sub(i),xlims);
		end
		if strcmp(ylims,'auto'),
			ylim(sub(i),'auto');
		elseif ~isempty(ylims),
			ylim(sub(i),ylims);
		end
	end
end
