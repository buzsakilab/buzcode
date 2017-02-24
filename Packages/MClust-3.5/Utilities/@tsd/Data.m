function v = Data(tsa, ix)

% ctsd/Data
%   d = Data(tsa)
%   d = Data(tsa, alignments)
%	Retrieves data from ctsd
%   if called with alignment list (timestamps), returns those tsa.Data(ix)
%   if called without, returns complete tsa.Data
%
% ADR 1998
% v4.1 30 oct 1998 correctly handles multi-d input
% version L4.1
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

switch nargin
case 2
   f = findAlignment(tsa, ix);
   v = SelectAlongFirstDimension(tsa.data,f);
  case 1
    v = tsa.data;
  otherwise
    error('Unknown number of input arguments');
end
