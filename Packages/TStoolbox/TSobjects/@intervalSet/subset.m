function A = subset(O, ix)

%  Returns subsets of an intervalSet
%  
%  	USAGE:
%  	A = subset(A, ix) 
%  	
%  	INPUTS:
%  	A  - an intervalSet object
%  	ix - index of the subsets
  
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  if nargin < 2 
    error('Call with two  arguments');
  end
  

  A = intervalSet(O.start(ix), O.stop(ix));
  
  