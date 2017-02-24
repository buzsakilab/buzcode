function tsa = ts(t,tUnits)
%
% TSD = ts(t,units)
% TSD = ts(t)
%
% A ts object contains a sequence of times. The units should be specified.
% If units are not specified, 'sec' are assumed.
%
% Methods
%    ts/Data         - Returns the timestamps as a matlab array
%    ts/StartTime    - First timestamp
%    ts/EndTime      - Last timestamp
%    ts/Restrict      - Keep timestamps within a certain range
%
% ADR
% version L5.0
% v5.0: JCJ 12/19/2002 includes support for time units 
% v5.1: JCJ 3/2/2003 allows convertion of tsd to ts. 
%
% Status: IN PROGRESS 

tsa.t = NaN;
tsa.units = NaN;

if nargin == 0
   tsa.t = [];
   tsa.units =[];
elseif nargin ==1
    if isa(t,'tsd')||isa(t,'ctsd')||isa(t,'ts')
        tsa.t=range(t);
        if isempty(strmatch('units',fieldnames(t)))
            warning('TS:units','t units not specified' )
        else
            tsa.units=units(t);
        end
        
    else
        tsa.t = t;
        tsa.units = 'sec';
    end
    
elseif nargin ==2
    if isa(t,'tsd')||isa(t,'ctsd')||isa(t,'ts')
        tsa.t=range(t,tUnits);
        tsa.units = tUnits;
    else
        tsa.t = t;
        tsa.units = tUnits;
    end
  
end

tsa = class(tsa, 'ts');
