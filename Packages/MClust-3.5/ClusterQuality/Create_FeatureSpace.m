function Fet = Create_FeatureSpace(WV)

% Fet = Create_FeatureSpace(WV)
%
% Creates an 8 dimensional feature space (energy and first princiapl
% component coefficient of the energy normalized waveform, calculated for 
% four tetrode channels) for use in calculating cluster quality values
%
% INPUTS: WV - a # of spikes X # of channels X # of waveform samples
%
% OUTPUTS: Fet - a # of spikes X # of features matrix
%                column 1 = energy on channel 1, column 2 = energy on
%                channel 2... column 5 = first principal component
%                coefficient on channel 1, etc. 

[nSpikes, nCh, nSamp] = size(WV);
Fet = repmat(NaN,nSpikes,8);

% Calculate energy
for iCh = 1:nCh
   Fet(:, iCh) = sqrt(sum(squeeze( WV(:, iCh, :) .* WV(:, iCh, :) ), 2));
end

% Normalize energy by the nubmer of waveform samples
Fet = Fet/nSamp;

% Calculate 1st principal component coefficients for the energy normalized
% waveforms 
for iCh = 1:nCh
    % pull out the waveforms from channel iCh
    w = squeeze(WV(:,iCh,:));    
    
    % normalize waveforms to unit L2 norm (so that only their SHAPE or
    % relative angles but not their length (energy) matters)
    l2norms = sqrt(sum(w.^2,2));
    w = w./l2norms(:,ones(1,nSamp));
    
    [pc,wpc,variances,t2] = princomp(w);
    
    % Get the first principal component coefficients for each spike
    Fet(:,iCh + 4) = wpc(:,1);
end