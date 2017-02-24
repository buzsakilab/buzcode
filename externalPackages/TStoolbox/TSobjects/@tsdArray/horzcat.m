function s = horzcat(varargin)

% Horizontal concatenation (overload of the [a, b] operator)  
% 
%  	USAGE:
%  	s = horzcat(tsa, tsb,..) 
%  
%	arguments must be all tsdArray.
%
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  
  C = {};
  
  for i = 1:length(varargin)
    if ~isa(varargin{i}, 'tsdArray')
      error('can only cat tsdArray');
    end
    A = varargin{i};
    C = [C A.C];
  end
  
  s = tsdArray(C);
  
  
