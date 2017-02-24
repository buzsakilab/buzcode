function C = compoundEvents(condTime, condSign, condSkip, condOffset)
% compoundEvents(condTime, condSign) constructor for the compoundEvents
% class
%  
% a compoundEvents object define a class of events formed by several
% events in a time series of events (e.g., single cell spike trains). 
% INPUTS
% condTime: a list of times representing limits the intervals between
% conditions must obey 
% condSign: defines how to interpret the limits -1 = less then, 0 = equal
% 1 = greater than
% condSkip: for each limit is n, and the limit will be applied to the
% interval between the i-th and the (i+n)-th event
% condOffset: if an event is recognized, the "center" of the event will
% be interpreted as the n-th spike in the event 
  
% copyright (c) 2004 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html

  if nargin == 0
    C.condTime = [];
    C.condSign = [];
    C.condSkip = [];
    C.condOffset = [];
    
    C = class(C, 'compoundEvents');
  end
  
  if(size(condTime, 2) ~= 1)
    error('condTime must be a column vector')
  end
  
  if(any(size(condSign) ~= size(condTime)))
    error('condSign must have the same size as condTime');
  end
  
  if any(condSign > 1 | condSign < -1)
    error('condSign may be only -1, 0, 1');
  end
  
  
  
  
  if(any(size(condSkip) ~= size(condTime)))
    error('condSkip must have the same size as condTime');
  end
  
  if any(condSkip < 1)
    error('condSkip may be only 1 or greater');
  end
  
  

  
  
  C.condTime = condTime;
  C.condSign = condSign;
  C.condSkip = condSkip;
  C.condOffset = condOffset;
  
  C = class(C, 'compoundEvents');
  
  