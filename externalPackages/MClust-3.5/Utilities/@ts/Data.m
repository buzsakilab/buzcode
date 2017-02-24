function t = Data(TS, ts)

% @ts/Data:
% t = Data(TS, ts)
% t = Data(TS)
%
% if ts included then returns the value just beneath ts, 
% if ts not included, returns an unclassed array of the data in TS

% ADR
% version L4.2
% v4.2, 18 nov 98 allowed ts flag.
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

if nargin == 1
   t = TS.t;
else
   t = TS.t(binsearch(TS.t, ts));
end
