function [ intcoupling ] = ISPCint(signal1,signal2,frange,int,sf)
%ISPCint(LFP,frange,int) Calculates the phase coherence (Inter-Site Phase
%Clustering) between two signals for a set of intervals
%
%
%%
%Filter the signals to get phase
[~,~,phase1] = FiltNPhase(signal1,frange,sf);
[~,~,phase2] = FiltNPhase(signal2,frange,sf);

%%
phase1 = IsolateEpochs2(phase1,int,0,sf);
phase2 = IsolateEpochs2(phase2,int,0,sf);

%%
numints = length(phase1);
intcoupling = zeros(numints,2);
for ii = 1:numints
    [intcoupling(ii,1),intcoupling(ii,2)] = ISPC(phase1{ii},phase2{ii});
end

