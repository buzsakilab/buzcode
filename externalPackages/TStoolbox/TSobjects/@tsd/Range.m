function R = Range(tsa, tsflag)

% Returns TSD timestamps
% 
%  USAGE:
%  R = Range(tsa)
%  R = Range(tsa, tsflag)
%  
%  OPTIONS:
%  tsflag - if 'ts' returns time in timestamps (default),
%  	    if 'sec' returns time in sec
%  	    if 'ms' returns time in ms


% ADR 
% version L4.1
% status: PROMOTED

% v4.1 28 oct 1998 flag no longer optional.
  if nargin == 2
    if isa(tsflag, 'units')
      ;
      
    elseif isa(tsflag, 'char')
      tsflag = time_units(tsflag);
    else
      error('tsflag must be time_units, or string');
    end
    
  else 
    tsflag = tsa.time_unit; %defaults to ts
  end
    
  cnvrt = convert(tsa.time_unit, tsflag);
  
  if cnvrt ~= 1
    R = tsa.t * cnvrt;
  else
    R = tsa.t;
  end

  