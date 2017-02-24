
function [redraw, rekey, undoable] =ShowAverageWaveform(iClust)

% [redraw, rekey, undoable] =ShowWaveformDensity(iClust)
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

redraw = false; rekey = false; undoable = false; % don't need to update

global MClust_TTData MClust_Clusters MClust_FeatureSources
global MClust_Colors
% MClust_AverageWaveform_ylim - added 3.5 ncst
global MClust_AverageWaveform_ylim

run_avg_waveform = 1;

[f MClust_Clusters{iClust}] = FindInCluster(MClust_Clusters{iClust});

% modified ncst 23 May 02
%         if length(f) > 150000
%             ButtonName=questdlg([ 'There are ' num2str(length(f)) ' points in this cluster. MClust may crash. Are you sure you want to check clusters?'], ...
%                 'MClust', 'Yes','No','No');
%             if strcmpi(ButtonName,'No')
%                 run_avg_waveform = 0;
%             end
%         else
if isempty(f)
    run_avg_waveform = 0;
    msgbox('No points in cluster.')
end

% modified to deal with waveforms that are not 32 samples and to
% use MClust_AverageWaveform_ylim
% - added 3.5 ncst
if run_avg_waveform
    clustTT = ExtractCluster(MClust_TTData, f);
    [mWV sWV] = AverageWaveform(clustTT);
    nWVSamples = size(mWV,2);
    AveWVFig = figure;
    for it = 1:4
        xrange = ((nWVSamples + 2) * (it-1)) + (1:nWVSamples);
        figure(AveWVFig);
        hold on;
        plot(xrange, mWV(it,:));
        errorbar(xrange,mWV(it,:),sWV(it,:)) %,'blue');
        set(findobj(gcf,'Type','line'),'Color',MClust_Colors(iClust + 1,:));
    end
    axis off
    axis([0 4*(nWVSamples+2) MClust_AverageWaveform_ylim])
    hold off
    title(['Average Waveform: Cluster ' num2str(iClust)]);
end
