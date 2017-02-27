function [d,dt]=formatts(d)
%
% Usage: [d,dt]=formatts(d)
%
% Helper function for CWT,XWT,WTC
% 
% Brings a timeseries into a shape so that it has two columns: [time, value].
% 
%
% (C) Aslak Grinsted 2002-2004
%

% -------------------------------------------------------------------------
%   Copyright (C) 2002-2004, Aslak Grinsted
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.



if (ndims(d)>2)
    error('Input time series should be 2 dimensional.');
end
if (numel(d)==length(d))
    d=[(1:length(d))' d(:)];
end
if size(d,1)<size(d,2)
    d=d';
end
if (size(d,2)~=2)
    error('Time series must have 2 columns.')
end

if (d(2,1)-d(1,1)<0)
    d=flipud(d);
end

dt=diff(d(:,1));
if any(abs(dt-dt(1))>1e-1*dt(1))
    error('Time step must be constant.');
end
if (dt==0)
    error('Time step must be greater than zero.')
end

dt=dt(1);

