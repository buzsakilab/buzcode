function TSO = find(C, TSA)
% TSO = find(C, TSA) find a compound event in a tsd, return tsd with
% event centers only.
%
% INPUTS:
% C: a compoundEvents object
% TSA: a tsd object
% OUTPUTS:
% TSO: a tsd object containing only the points relative to the event
% centers
  
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  
  % keyboard;
  ix = findIdx(C, TSA);
  
  TSO = subset(TSA, ix);
  
  