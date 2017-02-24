function [redraw, rekey, undoable] = ShowClusterInfo(iClust)
% [redraw, rekey,  undoable] = ShowClusterInfo(iClust)
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

global MClust_Clusters

msg{1} = GetName(MClust_Clusters{iClust});
msg{2} = ['class(' class(MClust_Clusters{iClust}) ')'];
msg{3} = '-----';
msg = cat(2,msg,GetInfo(MClust_Clusters{iClust}));
msgbox(msg, ['Cluster ', num2str(iClust)]);
