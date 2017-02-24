function [S,shank] = LoadMUAData(fbasename,varargin)

% Loads MUA data (cluster 1)
%
% USAGE:
% 	[S,shank] = LoadMUAData(fbasename)
%
% INPUTS:
%	fbasename: file base name
%	
% OUTPUTS:
%	S: tsdArray of MUA (one cell element per electrode group)
%	groups: electrode groups associated with each element
%
% OPTION:
%	[S,shank] = LoadMUAData(fbasename,groups)

% Adrien Peyrache, 2011


Fs = 20000;

if ~isempty(varargin)
    channels = varargin{1};
else
    syst = LoadXml([fbasename '.xml'],'raw');
    channels = (1:length(syst.SpkGrps));
end

try
    channels = channels(:)';
    S = {};
    shank = [];
    for ch=channels

        if exist([fbasename '.clu.' num2str(ch)],'file') && exist([fbasename '.res.' num2str(ch)],'file')
            clu = load([fbasename '.clu.' num2str(ch)]);

            clu = clu(2:end);
            if any(clu==1)
                res = load([fbasename '.res.' num2str(ch)]);
                S = [S;{ts(res(clu==1)/Fs)}];
                shank = [shank;ch];
            end
            
        end

    end

    S = tsdArray(S);
    
catch
    warning(lasterr);
    keyboard
end
