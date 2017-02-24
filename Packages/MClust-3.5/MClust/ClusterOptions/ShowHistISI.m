function [redraw, rekey, undoable] =ShowHistISI(iClust)
% [redraw, rekey, undoable] =ShowHistISI(iClust)
%
% INPUTS
%    iClust
%
% OUTPUTS
%
% NONE
% TO USE WITH MCLUST, put this in the MClust/ClusterOptions folder

% ADR 2003
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.
% Extensively modified by ADR to accomodate new ClusterOptions methodology
% v3.1 JCJ 2007: Added extra plot info
% ADR removed extraneous +1

redraw = false; rekey = false; undoable = false; % don't need to update

threshold=2; %2 ms

global MClust_Clusters MClust_FeatureSources MClust_FeatureTimestamps
[f MClust_Clusters{iClust}] = FindInCluster(MClust_Clusters{iClust});
[H,N]=HistISI(ts(MClust_FeatureTimestamps(f)),'DoPlotYN','yes');
a=axis;
hold on
hold on;plot([1 1]*threshold,[a(3:4)],'g:')       
title(['Cluster ' num2str(iClust) ' ISI Histogram:  ' num2str(sum(H(N<threshold))) '/' num2str(sum(H)) ' ISIs less than ' num2str(threshold) 'ms']);
        
