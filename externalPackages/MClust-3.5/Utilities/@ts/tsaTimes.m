function S = tsaTimes(tsa,unitflag,opflag)
% ctsd/Range
% 
%  R = tsaTimes(tsa,unitflag,opflag)
%
%  returns times requested by opflag in units of unitflag
%
% INPUTS
%     tsa         
%     unitflag:  if 'ts' returns time in timestamps,
%                if 'sec' returns time in sec
%                if 'sec0' returns time in sec counting from 0
%                if 'ms' returns time in ms
%
%     opflag:    if 'a' returns all timestamps converted
%                if 's' returns starting timestamp (1st)
%                if 'e' returns ending timestamp (last)
%                if 'd' returns DT (difference between timestamps)
%
% JCJ
% version 1.0


switch (opflag)
    case 'a'
        S = tsa.t;
        
    case 's'
        S = min(tsa.t);
        
    case 'e'
        S = max(tsa.t);
        
    case 'd'
        S = nanmedian(diff(tsa.t));

    otherwise
        error('Invalid opflag: no conversion possible');
end
        

if isempty(strmatch('units',fieldnames(tsa)))
    warning('TS:UNITS','tsa units not specified: assuming tsa units = sec ' )
    S = convert(S,'sec',unitflag);
    return
else
    S = convert(S,tsa.units,unitflag);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = convert(R,unit,unitflag)

switch (unitflag)
    case 'sec'
        switch unit
            case 'sec'
                % R = R;
            case 'ts'
                R = R/10000;
            case 'ms'
                R = R/1000;
            otherwise
                warning('TS:Units','tsa has invalid units: no conversion possible' );
        end
        
    case 'sec0'
        switch unit
            case 'sec'
                R = (R - min(R));
            case 'ts'
                R = (R - min(R))/10000;
            case 'ms'
                R = (R - min(R))/1000;
            otherwise
                warning('TS:units','tsa has invalid units: no conversion possible' );
        end
        
    case 'ts'
        switch unit
            case 'sec'
                R = R*10000;
            case 'ts'
                % R = R;
            case 'ms'
                R = R*10;
            otherwise
                warning('TS:units', 'tsa has invalid units: no conversion possible' );
        end
        
    case 'ms'
        switch unit
            case 'sec'
                R = R*1000;
            case 'ts'
                R = R/10;
            case 'ms'
                % R = R;
            otherwise
                warning('TS:units', 'tsa has invalid units: no conversion possible' );
        end
    otherwise
        warning('TS:units', 'Convert called with invalid unitflag: no conversion possible');
end
