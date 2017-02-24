function [S,shank,cellIx] = LoadSpikeData(fbasename,varargin)

% Load spikes times from a clu and res file
%
% USAGE:
%     [S,groups,cellIx] = LoadSpikeData(fbasename)
%
% INPUT:
%    fbasename: session file basename
%
% OUTPUT:
%    S: a tsdArray of spike trains
%    groups: a vector giving electrode group for each neuron
%    cellIx: a vector giving cell's clu index in the clu file

% Adrien Peyrache, 2011
syst = LoadXml([fbasename '.xml']);
Fs = syst.SampleRate; 

if ~isempty(varargin)
    channels = varargin{1};
else
    channels = (1:length(syst.SpkGrps));
end

try
    channels = channels(:)';
    S = {};
    shank = [];
    cellIx = [];
    for ch=channels

        if exist([fbasename '.clu.' num2str(ch)],'file') && exist([fbasename '.res.' num2str(ch)],'file')
            clu = load([fbasename '.clu.' num2str(ch)]);
            res = load([fbasename '.res.' num2str(ch)]);

            clu = clu(2:end);
            cluIx = unique(clu);
            cluIx = cluIx(cluIx>1);
            cluIx = cluIx(:)';
                for c=cluIx
                    S = [S;{ts(res(clu==c)/Fs)}];
                    shank = [shank;ch];
                    cellIx = [cellIx;c];
                end
        end

    end

    S = tsdArray(S);
    
catch
    warning(lasterr);
    keyboard
end
