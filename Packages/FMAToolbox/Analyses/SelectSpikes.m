function selected = SelectSpikes(spikes,varargin)

%SelectSpikes - Discriminate bursts vs single spikes.
%
%  USAGE
%
%    selected = SelectSpikes(spikes,<options>)
%
%    spikes         spike times
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'mode'        either 'bursts' (default) or 'single'
%     'isi'         max inter-spike interval for bursts (default = 0.006),
%                   or min for single spikes (default = 0.020)
%    =========================================================================
%
%  OUTPUT
%
%    selected       a logical vector indicating for each spike whether it
%                   matches the criterion

% Copyright (C) 2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
isiBursts = 0.006;
isiSpikes = 0.020;
mode = 'bursts';

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help SelectSpikes">SelectSpikes</a>'' for details).');
end

% Check parameter size
if ~isdvector(spikes),
	error('Incorrect spikes (type ''help <a href="matlab:help SelectSpikes">SelectSpikes</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
  if ~ischar(varargin{i}),
    error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help SelectSpikes">SelectSpikes</a>'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'mode',
      mode = varargin{i+1};
      if ~isstring_FMAT(mode,'bursts','single'),
        error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help SelectSpikes">SelectSpikes</a>'' for details).');
      end
    case 'isi',
      isi = varargin{i+1};
      if ~isdscalar(isi,'>0'),
        error('Incorrect value for property ''isi'' (type ''help <a href="matlab:help SelectSpikes">SelectSpikes</a>'' for details).');
      end
    otherwise,
      error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SelectSpikes">SelectSpikes</a>'' for details).']);
  end
end

t = spikes(:);
dt = diff(t);

if strcmp(mode,'bursts'),
	b = dt<isi;
	selected = logical([0;b])|logical([b;0]); % either next or previous isi < threshold
else
	s = dt>isi;
	selected = logical([0;s])&logical([s;0]); % either next or previous isi < threshold
end
