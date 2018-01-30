function [ phaseamplitudemap,ampfreqs,phasecenters ] = PhaseAmplitudeDist( LFP,phaserange,amprange,int,sf_LFP )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% DEV
%datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/';
%recname = '20140526_277um';
% recname = 'Dino_061814_mPFC';
 %figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/TransitionsAndSpindles/AnalysisFigsMats/Spindle LFP/';
% 
 %load([datasetfolder,recname,'/',recname,'_LFP.mat'])
 %load([datasetfolder,recname,'/',recname,'_StateIntervals.mat'])
% 
% 
 %int = StateIntervals.Spindles;
 %sf_LFP = 1250;
% 
 %phaserange = [10 20];
 %amprange = [24 256];

nfreqs = 100;
ncyc = 7;

%% Deal with input types

if isa(int,'intervalSet')
    int = [Start(int,'s'), End(int,'s')];
end

t_LFP = (1:length(LFP))'/sf_LFP;
%LFP = ZScoreToInt(LFP,int);
%% Filter LFP
[~,~,LFP_phase] = FiltNPhase(LFP,phaserange,sf_LFP);


%% Get LFP, Phase in intervals
edgebuffer = 1; %s
edgebuffer_si = edgebuffer.*sf_LFP;
edgebuffer = edgebuffer.*[1 1];

LFP_int = IsolateEpochs2(LFP,int,edgebuffer,sf_LFP);
LFPphase_int = IsolateEpochs2(LFP_phase,int,edgebuffer,sf_LFP);

%% Wavelet Transform LFP in intervals
[ampfreqs,~,spec_int] = WaveSpec(LFP_int,amprange,nfreqs,ncyc,1/sf_LFP,'log');
spec_int = cellfun(@(X) abs(X),spec_int,'UniformOutput',false);
%Remove Buffers
spec_int = cellfun(@(X) X(:,edgebuffer_si:end-edgebuffer_si),spec_int,'UniformOutput',false);
LFPphase_int = cellfun(@(X) X(edgebuffer_si:end-edgebuffer_si),LFPphase_int,'UniformOutput',false);

%% Bin phase and power
numbins = 100;
phasebins = linspace(-pi,pi,numbins+1);
phasecenters = phasebins(1:end-1)+(phasebins(2)-phasebins(1));

specall = cat(2,spec_int{:})';
phaseall = cat(1,LFPphase_int{:});
[phasedist,~,phaseall] = histcounts(phaseall,phasebins);


phaseamplitudemap = zeros(numbins,nfreqs);
for bb = 1:numbins
    phaseamplitudemap(bb,:) = mean(specall(phaseall==bb,:),1);
end
    
%% Plot
figure
    imagesc(phasecenters,log2(ampfreqs),zscore(phaseamplitudemap)')
    LogScale('y',2)
    axis xy

end

