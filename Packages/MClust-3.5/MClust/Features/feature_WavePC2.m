function [wavePCData, wavePCNames, wavePCPar] = feature_WavePC2(V, ttChannelValidity, Params)

% MClust
% [Data, Names] = feature_wavePCA(V, ttChannelValidity)
% Calculate first waveform PRINCIPAL COMPONENTS  (PC2)
%
% INPUTS
%    V = TT tsd
%    ttChannelValidity = nCh x 1 of booleans
%    Params = 4x1 CellArray struct with fields
%             Params{}.pc (eigenvectos),
%             Params{}.av (averages),
%             Params{}.sd (std deviations)
%
% OUTPUTS
%    Data - nSpikes x nPC*nCh  of waveform PCs of each spike for each valid channel
%    Names - "wavePCn: Ch"

% CCC: PL
% version 1.0
%
% Status: PROMOTED (Release version)
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

%%% PARAMETERS:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iPC  = 2;    % number of principal component to keep (per channel)
norm = 1;    % normalize Waveforms (1) or don't normalize (0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 2
        recalcParams = 1;
    case 3
        recalcParams = 0;
    otherwise
        error('feature_wavePC2 needs 2 or 3 input arguments');
end


TTData = Data(V);
[nSpikes, nCh, nSamp] = size(TTData);

f = find(ttChannelValidity);
lf = length(f);

wavePCNames = cell(lf, 1);
wavePCData  = zeros(nSpikes, lf);
I = ones(nSpikes,1);

for iC = 1:lf
    w = squeeze(TTData(:,f(iC),:));    % get data in nSpikes x nSamp array
    if norm
        % normalize waveforms to unit L2 norm (so that only their SHAPE or
        % relative angles but not their length (energy) matters)
        l2norms = sqrt(sum(w.^2,2));
        w = w./l2norms(:,ones(1,nSamp));
        % ncst 3.5 modification - if you have a waveform of all zeros
        % (encountered in some Axona recordings on valid channels), you'll get
        % a NaN after normalization, which kills the princomp calculation.  So,
        % replace NaNs with zeros.
        w(isnan(w(:,1)),:) = 0;
    end
    if recalcParams
        cv = cov(w);
        sd = sqrt(diag(cv))';        % row std vector
        av = mean(w);                % row mean vector
        pc = wavePCA(cv);            % get PCA eigenvectors (in columns of pc)
    else
        pc = Params{f(iC)}.pc;
        av = Params{f(iC)}.av;
        sd = Params{f(iC)}.sd;
    end
    wavePCPar{f(iC)}.pc = pc;
    wavePCPar{f(iC)}.av = av;
    wavePCPar{f(iC)}.sd = sd;
    wstd=(w-(I*av))./(I*sd);     % standardize data to zero mean and unit variance
    wpc = wstd*pc;               % project data onto principal component axes
    wavePCData(:,iC) = wpc(:,iPC);
    wavePCNames{iC} = ['WavePC2'  ': ' num2str(f(iC))];
end
