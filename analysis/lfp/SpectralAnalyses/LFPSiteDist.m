function [ sitesim ] = LFPSiteDist(cohmat,spikegroups,probegroups)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
numprobes = length(probegroups);
allsites = [spikegroups{:}];

onprobesim = zeros(size(cohmat));
offprobesim = zeros(size(cohmat));
sitesim = zeros(size(cohmat));
for pp = 1:numprobes
    probesites = spikegroups(probegroups{pp});
    probesites = [probesites{:}];    
    othersites = setdiff(allsites,probesites);
    
    onprobesim(probesites+1,probesites+1) = 1-squareform(pdist(cohmat(probesites+1,probesites+1),'cityblock'))/length(probesites);

    offprobesim(probesites+1,probesites+1) = 1-squareform(pdist(cohmat(probesites+1,othersites+1),'cityblock'))/length(probesites);
    
    sitesim(probesites+1,probesites+1) = cohmat(probesites+1,probesites+1)+onprobesim(probesites+1,probesites+1)+offprobesim(probesites+1,probesites+1);
end


%% Figure
% figure
%     subplot(2,3,1)
%         imagesc(cohmat([spikegroups{:}]+1,[spikegroups{:}]+1))
%         colorbar
%     subplot(2,3,2)
%         imagesc(onprobesim([spikegroups{:}]+1,[spikegroups{:}]+1))
%         colorbar
%         
%     subplot(2,3,3)
%         imagesc(offprobesim([spikegroups{:}]+1,[spikegroups{:}]+1))
%         colorbar
%         
%     subplot(2,3,4)
%         imagesc(sitesim([spikegroups{:}]+1,[spikegroups{:}]+1)./3)
%         colorbar
end

