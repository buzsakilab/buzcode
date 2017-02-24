function sortedCA = sortcell(stringCA, k)

% sortedCA = sortcell(stringCA)
%
% sorts a cell array of strings
%
% INPUTS
%    stringCA = a cell array of strings
%    k (opt) = line[s] on which to sort on
%
% OUTPUTS
%    sortedCA = a cell array of strings
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

if nargin == 1
   sortedCA = cellstr(sortrows(char(stringCA)));
else
   sortedCA = cellstr(sortrows(char(stringCA), k));
end