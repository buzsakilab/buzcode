function [redraw, rekey, undoable] = EvalOverlap

% [redraw, rekey, undoable] = EvalOverlap
%
% INPUTS
%
% OUTPUTS
%
% NONE
% TO USE WITH MCLUST, put this in the MClust/ClusterOptions folder
%
% ADR 2008
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.5.
% Version control M3.5.

redraw = false; rekey = false; undoable = false;

    global MClust_Clusters MClust_FeatureTimestamps
    nS = length(MClust_FeatureTimestamps);
    nClust = length(MClust_Clusters);
    
    % ncst added 16 may 02
    if nClust == 1
        errordlg('Only one cluster exists.', 'MClust error', 'modal');
        return
    end
    
    nToDo = nClust * (nClust-1)/2;
    iDone = 0;
    overlap = zeros(nClust);
    for iC = 1:nClust
        [fI] = FindInCluster(MClust_Clusters{iC});
        iSpikes = false(nS,1); iSpikes(fI) = true;
        for jC = (iC+1):nClust
            iDone = iDone +1;
            DisplayProgress(iDone, nToDo, 'Title', 'Evaluating overlap');
            [fJ] = FindInCluster(MClust_Clusters{jC});
            jSpikes = false(nS,1); jSpikes(fJ) = true;
            overlap(iC,jC) = sum(iSpikes & jSpikes);
            overlap(jC,iC) = overlap(iC,jC);
        end
    end
    
    overlapT = [(0:nClust)', [1:nClust; overlap]];
        
    figure; h = bar3(overlap); title('Overlap');
    
    msgbox(num2str(overlapT), 'Overlap');
    
