function SNR = GetSNR(iClust,varargin)

% GetSNR(iClust)
%
% INPUTS
%    iClust
%
% OUTPUTS
%
% SNR - signal to noise ratio of each channel for the cluster
%       Calculated for each channel as the value of the peak (or valley if
%       the cell is inverted) minus the value of the first waveform sample
%       divided by the standard deviation of the first sample.  Calculation
%       based on a poster ncst saw at SFN 2003.  
%
%  If no outputs are requested, displays the maximum SNR of all channels
%
% TO USE WITH MCLUST, put this in the MClust/ClusterOptions folder

% ADR 2003
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.
% Extensively modified by ADR to accomodate new ClusterOptions methodology

global MClust_Clusters MClust_FeatureData MClust_TTdn MClust_TTfn MClust_ChannelValidity MClust_TTData MClust_FeatureTimestamps

ClustTT = [];
NoiseTT = [];

Extract_varargin;

if isempty(ClustTT) || isempty(NoiseTT)

	[f MClust_Clusters{iClust}] = FindInCluster(MClust_Clusters{iClust}, MClust_FeatureData);
	
	ClustTT = ExtractCluster(MClust_TTData,f);
	NoiseTT = ExtractCluster(MClust_TTData,1:length(MClust_FeatureTimestamps));
end

mWV = AverageWaveform(ClustTT); 
[Noise_mWV Noise_sWV] = AverageWaveform(NoiseTT); 
[V maxPeak] = max(abs(mWV)');

Peak_vals = mWV(sub2ind(size(mWV),1:4,maxPeak))';
Peak_value = abs(Peak_vals) - sign(Peak_vals).*Noise_mWV(:,1);
Noise_std = Noise_sWV(:,1);

SNR = Peak_value./Noise_std;

if nargout == 0
	disp(sprintf('  Cluster %2.0f Signal-To-Noise ratio = %5.1f',iClust,max(SNR)));
end