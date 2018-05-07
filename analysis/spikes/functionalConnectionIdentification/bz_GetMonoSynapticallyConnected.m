function mono_res = bz_GetMonoSynapticallyConnected(basepath,varargin)

%% Method for estimating likelihood of monosynaptic connectivity based off of relative spike timing
% Based off of Stark and Abeles 2009 and verified for PYR -INT in English et al., 2017
% Assumes network synchrony is low frequency and that synapses induce high
% frequency additive "injected" excess synchrony above that low frequency
% network co-modulation. In addition this model assumes that synaptic
% contribution should be delayed 0.8 -2.8ms (this will change according to
% cell type and inter somatic distance!). 

% to get probabilities derived from English et al., 2017, you must have 'ProbSynMat.mat'


% see bz_MonoSynConvClick 

%%%  OPTIONAL INPUTS:
%%%
%%%  binSize = timebin to compute CCG (in seconds)
%%%
%%%  duration = window to compute CCG (in seconds)
%%%
%%%  epoch = [start end] (in seconds)
%%%
%%%  cells = N x 2 matrix -  [sh celID] to include (NOTE indexing will be
%%%          done on full spikeIDlist
%%%
%%%  conv_w = # of time bins for computing the CI for the CCG.
%%%
%%%  alpha = type I error for significance testing using convolution
%%%     technique. Stark et al, 2009
%%%
%%%  calls: CCG, InInterval,FindInInterval (from FMA toolbox)
%%%         tight_subplot, mtit (from matlabcentral)
%%%         bz_cch_conv, bz_PlotMonoSyn, bz_MonoSynConvClick
%%%
%%%  OUTPUT
%%%  mono_res.alpha = p-value
%%%  mono_res.ccgR = 3D CCG (time x ref x target;
%%%  mono_res.sig_con = list of significant CCG;
%%%  mono_res.Pred = predicted Poisson rate;
%%%  mono_res.Bounds = conf. intervals of Poisson rate;
%%%  mono_res.conv_w = convolution windows (ms)
%%%  mono_res.completeIndex = cell ID index;
%%%  mono_res.binSize = binSize;
%%%  mono_res.duration = duration;
%%%  mono_res.manualEdit = visual confirmation of connections
%%%  mono_res.Pcausal = probability of getting more excess in the causal than anticausal direction;
%%%  mono_res.FalsePositive = FalsePositive rate from English et al., 2017;
%%%  mono_res.TruePositive = TruePositive rate from English et al., 2017;

%Written by Sam McKenzie 2017

%load data
spikes = bz_GetSpikes('basepath',basepath);

%get shank clu
spikeIDs = [spikes.shankID(spikes.spindices(:,2))' spikes.cluID(spikes.spindices(:,2))' spikes.spindices(:,2)];

%call main script
mono_res = bz_MonoSynConvClick (spikeIDs,spikes.spindices(:,1),varargin);
end