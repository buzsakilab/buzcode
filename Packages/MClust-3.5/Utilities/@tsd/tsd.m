function tsa = tsd(t, qData, tUnits)
%
% tsa = tsd(t,data)
% tsa = tsd(t,data, units)  e.g. units = {'ts','sec', or 'ms'; assumes 'sec'}
%
% tsd is a class of "timestamped arrays"
% 	It includes a list of timestamps
% 	and data (possibly an array).  
%    The first dimension of Data correspond to the
%    timestamps given.  The timestamps must be sequential,
% 	but do not have to be continuous.
%
% Methods
%    tsd/Range     - Timestamps used
%    tsd/Data      - Returns the data component
%    tsd/DT        - Returns the DT value (mean diff(timestamps))
%    tsd/StartTime - First timestamp
%    tsd/EndTime   - Last timestamp
%    tsd/Restrict  - Keep data within a certain range
%    tsd/Mask      - Make all non-mask values NaN
%    tsd/AntiMask  - Make all mask values NaN
%
%  Note: data can be N-dimensional, time is always the first axis.
%
% ADR 
% version L6.0
%  
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0+.
% Version control M3.5.

switch nargin
   
case 0
 tsa.t = [];
 tsa.data = [];
 tsa.units= [];


case 1
   if isa(t, 'tsd')
      tsa = t;
      return;
   elseif isa(t, 'ctsd')
      tsa.t = Range(t);
      tsa.data = Data(t);
      if ~isfield(struct(t), 'units')
          warning('TSD:Units','units not specified' )
      else
          tsa.units=Units(t);
      end      
  elseif isa(t,'struct')
      tsa.t=t.t;
      tsa.data=t.data;
      if ~isfield(t, 'units')
          warning('TSD:Units','units not specified' )
      else
          tsa.units=t.units;
      end
      
  else
      error('Unknown copy-from object');
  end
   
case 2
   tsa.t = t;
   tsa.data = qData;
   tsa.units= 'sec';
case 3
   tsa.t = t;
   tsa.data = qData;
   tsa.units= tUnits;   
otherwise
   error('Constructor error tsd');
end

tsa = class(tsa, 'tsd');
   
