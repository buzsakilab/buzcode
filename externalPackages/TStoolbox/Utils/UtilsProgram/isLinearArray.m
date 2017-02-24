function t = isLinearArray(X)
% t = function isLinearArray(X) returns one if X is row or column vector
% 
% 

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html



t = 0

if ~isa(X, 'numeric')
  return
end

if length(X) > 2
  return
end

t = any(size(X) == 1);

 