function tss = subset(tsa, ix);

%  Returns tsd of a subset of input
%  
%  	USAGE:
%  	tss = subset(tsa, ix); 
%  	
%  	INPUTS:
%  	ts - an input object
%  	ix - an array of data indices 
%  	
%  	OUTPUTS: 
%  	tss - a tsd containing the columns indicated by ix

% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
  
 
  t = Range(tsa);
  d = Data(tsa);
  if max(ix) > size(d,2)
      error('Out of range')
  end
  tss = tsd(t, d(:,ix));
  
  