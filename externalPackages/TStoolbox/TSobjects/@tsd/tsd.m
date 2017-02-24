function tsa = tsd(t, Data, varargin)
  
% Constructor of TSD object
%  
%  	USAGE:
%  	tsa = tsd(t,data, 'optionName1', optionValue1,... )
%  
%  tsd is a class of "timestamped arrays"
%  composed of an array of timestamp (column vector), and an array of data
%  associated to each timestamp (first index corresponding to time).
%  data maybe empty, for 0-dimensional data (for example ts
%  data). Othewise first dimension must be same as length of timestamp
%  vector
%  timestamps are assumed to be sorted. If this is not the case, the
%  constructor failes, unless the fixOrder option is set
%  options may include the following:
%  'TimeUnits' - the unit in which the time is expressed (defaults to
%  		 ts=1/10000s)
%  'DataUnits' - the unit in which data are expressed (for the moment same
%  		 for all the data) 
%  'fixOrder'	- if set to nonzero, constructor will reorder the timestamps,
%  		  in the case they are not sorted.

% ADR 
% version L4.0
% status: PROMOTED

% version by Francesco P. Battaglia 2004  
  

if nargin == 0
 tsa.time_unit = time_units('ts');
 tsa.data_unit = units('none','none');
 tsa.t = NaN;
 tsa.data = NaN;
 tsa.name = '';
 tsa = class(tsa, 'tsd');
 return
end 


if nargin < 2
  error('tsd constructor must be called with T, Data');
end


% default option values
fixOrder = 0;

if ~isempty(varargin)
    defined_options = ...
    dictArray({ { 'TimeUnits', { time_units('ts'), {'units', 'char' } } }, 
		{ 'DataUnits', { units('none', 'none'), { 'units' } } },
		{ 'fixOrder', { 0, { 'numeric' } } },
        { 'keepUnordered', { 0, { 'numeric' } } },
	        { 'Name', { '', { 'char' } } } } );     

    opt_varargin = varargin;
    getOpt;
else
    TimeUnits = time_units('ts');
    DataUnits = [];
    fixOrder = 0;
    keepUnordered = 0;
    Name = '';
end
		  
tsa.time_unit = time_units(TimeUnits);
tsa.data_unit = DataUnits;


if isempty(t)
  t = zeros(0,1);
end

  
if all(size(t) ~= 1) 
  error('t must b row or column error');
end

if size(Data, 1) ~= length(t) & (~isempty(Data))
  error('Data must be N x X x Y x ..., where N is length(t), or empty');
end


t = t(:);

if any(diff(t) < 0)
  if fixOrder
    [t, ix] = sort(t);
    Data = SelectAlongFirstDimension(Data, ix);
  elseif ~keepUnordered
    error('timestamps are not sorted!');
  end
end




tsa.t = t;
tsa.data = Data;
tsa.name = Name;
tsa = class(tsa, 'tsd');

