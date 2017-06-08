function filtered = bz_Filter(samples,varargin)

%Filter - Filter samples.
%
%  USAGE
%
%    filtered = bz_Filter(samples,<options>)
%
%INPUT
%    samples        samples given as a [t x channels] timeseries vector
%                   or as buzcode structure with fields samples.data,
%                   samples.timestamps, samples.samplingRate
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'passband'    pass frequency range
%     'stopband'    stop frequency range
%     'order'       filter order (number of cycles, default = 4)
%     'ripple'      filter ripple (default = 20)
%     'filter'      choose filter type between 'cheby2' (default) and 'fir1'
%     'nyquist'     nyquist frequency (default = 625), calculated 
%                   automatically from samples.samplingRate for BUZCODE
%     'FMAlegacy'   true/false, uses FMA legacy input style
%                   (samples given as a list of (t,v1,v2,v3...) tuples)
%    =========================================================================
%
%OUTPUT
%   filtered        -if the input is a timeseries vector, output is as well
%                   -if the input is a buzcode structure, output is a
%                    buzcode structure with the following fields:
%       .data       the filtered data
%       .phase      phase, calculated by the hilbert method
%       .amp        amplitude, calculated by the hilbert method
%       .timestamps     (s)
%       .samplingRate   (Hz)
%       .filterparms    a structure that holds the parameters used for
%                       filtering, for future reference
%
% Copyright (C) 2004-2011 by Michaël Zugaro
% updated 2017 DLevenstein for buzcode
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%% Parms and Defaults 

% Default values
passband = [];
stopband = [];
order = 4;
ripple = 20;
nyquist = 625; %written over later from samples.samplingRate if available
type = 'cheby2';
FMAlegacy = false;
BUZCODE = false;

% Check number of parameters
if nargin < 1 | mod(length(varargin),2) ~= 0,
	error('Incorrect number of parameters (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
end

if isempty(samples),
	filtered = samples;
	return;
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i) ' is not a property (type ''help <a href="matlab:help Filter">Filter</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'passband',
			if ~isempty(stopband),
				error('Cannot specify both a passband and stopband (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
			end
			passband = varargin{i+1};
			if ~isdvector(passband,'#2','>=0'),
				error('Incorrect value for ''passband'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
			end

		case 'stopband',
			if ~isempty(passband),
				error('Cannot specify both a passband and stopband (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
			end
			stopband = varargin{i+1};
			if ~isdvector(stopband,'#2','>=0'),
				error('Incorrect value for ''stopband'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
			end

		case 'filter',
			type = lower(varargin{i+1});
			if ~isstring_FMAT(type,'cheby2','fir1'),
				error(['Unknown filter type ''' type ''' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).']);
			end

		case 'order',
			order = lower(varargin{i+1});
			if ~isiscalar(order,'>0'),
				error('Incorrect value for ''order'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
			end

		case 'ripple',
			ripple = lower(varargin{i+1});
			if ~isiscalar(ripple,'>0'),
				error('Incorrect value for ''ripple'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
			end

		case 'nyquist',
			nyquist = varargin{i+1};
			if ~isiscalar(nyquist,'>0'),
				error('Incorrect value for property ''nyquist'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
            end
            
        case 'fmalegacy'
            FMAlegacy = varargin{i+1};
            if ~islogical(FMAlegacy)
                error('Incorrect value for property ''FMALegacy''');
            end

		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).']);

	end
end
%%
%Check the input structure for buzcode stuff
if isstruct(samples)
    BUZCODE = true;
    if ~isfield(samples,'data')
        error('Your input structure has no .data field :(')
    end
    if ~isfield(samples,'samplingRate')
        warning('Your input strucutre has no .samplingRate field, using provided/default nyquist frequency')
        samples.samplingRate = 2.*nyquist;
    end
    if ~isfield(samples,'timestamps')
        warning('Your input structure has no .timestamps field, one will be provided for you')
    end
    nyquist = 0.5.*samples.samplingRate;
end

if isempty(passband) && isempty(stopband),
	error('Missing passband or stopband (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
end


switch(type),
	case 'cheby2',
		if ~isempty(passband),
			if passband(1) == 0,
				[b a] = cheby2(order,ripple,passband(2)/nyquist,'low');
			else
				[b a] = cheby2(order,ripple,passband/nyquist);
			end
		else
			[b a] = cheby2(order,ripple,stopband/nyquist,'stop');
		end
	case 'fir1',
        %Order input to fir1 needs to be in samples, not cycles
		if ~isempty(passband),
			if passband(1) == 0,
                filt_order = round(order*2*nyquist./passband(2));    
				[b a] = fir1(filt_order,passband(2)/nyquist,'low');
            else
                filt_order = round(order*2*nyquist./passband(1));  
				[b a] = fir1(filt_order,passband/nyquist);
			end
        else
            filt_order = round(order*2*nyquist./stopband(1));
			[b a] = fir1(filt_order,stopband/nyquist,'stop');
		end
end


if FMAlegacy %FMA has (samples given as a list of (t,v1,v2,v3...) tuples)
    filtered(:,1) = samples(:,1);
    for i = 2:size(samples,2),
        filtered(:,i) = FiltFiltM(b,a,double(samples(:,i)));
    end
elseif BUZCODE %BUZCODE has samples as a data structure
    filtered.timestamps = samples.timestamps;
    for i = 1:size(samples.data,2),
        filtered.data(:,i) = FiltFiltM(b,a,double(samples.data(:,i)));
        filtered.amp(:,i) = abs(hilbert(filtered.data(:,i)));
        filtered.phase(:,i) = angle(hilbert(filtered.data(:,i)));
    end
    filtered.filterparms.passband = passband;
    filtered.filterparms.stopband = stopband;
    filtered.filterparms.order = order;
    filtered.samplingRate = samples.samplingRate;
else %or if you just want filter a basic timeseries
    for i = 1:size(samples,2),
        filtered(:,i) = FiltFiltM(b,a,double(samples(:,i)));
    end
end

end