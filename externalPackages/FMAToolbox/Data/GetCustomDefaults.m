function output = GetCustomDefaults(current,property,value)

%GetCustomDefaults - Get custom default value for a given function property.
%
%  This is an 'internal' function used by FMAToolbox. You should not need
%  to use it, unless you are developping new functions for this toolbox.
%
%  See Example below.
%
%  WARNING
%
%    Custom defaults should be allowed with great care, because they *implicitly*
%    change the default behavior of functions: because the actual function call
%    does not make it explicit that default values are being overriden, the same
%    code can yield different results for different users, or even for the same
%    user before and after she/he changes the custom defaults.
%
%  USAGE
%
%    output = GetCustomDefaults(current,property,value)
%
%    current        current value
%    property       property name
%    value          default value if custom default value is undefined
%
%  EXAMPLE
%
%    The function GetPositions has optional property-value pairs, including:
%
%     - 'mode' (default value 'clean')
%     - 'coordinates' (default value 'normalized')
%     - 'pixel' (no default value)
%     - 'distances' (default value [0 Inf])
%
%    Users may want to change default values, e.g. have 'coordinates'
%    set to 'video' whenever they call GetPositions without an explicit
%    value. GetPositions therefore includes the following code:
%
%      % Default values (customizable defaults must be empty at this point)
%      c = '';
%
%      % Parse parameters
%      (possibly set a value for c, if supplied among optional input parameters)
%
%      % Customizable defaults
%      c = GetCustomDefaults(c,'coordinates','normalized');
%
%    GetCustomDefaults will first test if c was already set (non-empty), in which
%    case c will remain unaltered. Otherwise, it will test if there is a user-defined
%    default value for 'coordinates', and issue a warning and set c to this value if
%    it exists, or set it to 'normalized' otherwise.
%
%    To define a custom default value, users would have to e.g. add the following
%    lines to their startup.m file:
%
%     global SETTINGS;
%     SETTINGS.GetPositions.coordinates = 'video';
%
%  NOTE
%
%    Do not forget to document when the default value for a property can be customized.
%    See the code for GetPositions for an example.


% Copyright (C) 2010-2014 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if nargin == 2,
	% Backward compatibility for former usage: output = GetCustomDefaults(property,value)
	value = property;
	property = current;
	warning('Deprecated usage (type ''help <a href="matlab:help GetCustomDefaults">GetCustomDefaults</a>'' for details).');
else
	% Non-empty input
	if ~isempty(current),
		output = current;
		return
	end
end

% Default value
output = value;

% Find out calling function name
stack = dbstack;
if length(stack) < 2, error('Stack error (sorry, cannot be more explicit)'); end
functionName = stack(2).name;

% Look for custom default
global SETTINGS;
if ~exist('SETTINGS'), return; end
if ~isfield(SETTINGS,functionName), return; end
p = getfield(SETTINGS,functionName);
if ~isfield(p,property), return; end
output = getfield(p,property);

% Warn user
if isdvector(output),
	d =  sprintf('%g ',output);
	d = d(1:end-1);
	if length(output) ~= 1,
		d = ['[' d '] '];
	else
		d = [d ' '];
	end
elseif isstring_FMAT(output),
	d = ['''' output ''' '];
else
	d = '';
end
warning(['Using custom default value ' d 'for ''' property ''' in function ''' functionName '''.']);

