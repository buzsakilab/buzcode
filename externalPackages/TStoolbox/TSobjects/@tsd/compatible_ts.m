function t = compatible_ts(tsa, tsb, tl)

% Tells whether two tsd's have the same timestamps, within relative tolerance
%  
%  	USAGE:
%  	t = compatible_ts(tsa, tsb, tol) 
%  	
%  	INPUTS:
%  	tsa,tsb: two tsd to compare
%  	tol: tolerance (default 0)
%  
%  	OUPUTS:
%  	t: scalar value
%  
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  
  
  tol = 0;
  
  if nargin > 2
    tol = tl;
  end
  
  
  t = 0;
  
  if length(tsa.t) ~= length(tsb.t)
    return
  end
  
  t = all((2 * abs(tsa.t -tsb.t) ./ (abs(tsa.t) + abs(tsb.t)+eps)) < tol );
  
  