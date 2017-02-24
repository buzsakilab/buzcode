function [redraw, rekey, undoable] = ShowXCorr(iClust)

% [redraw, rekey, undoable] =ShowXCorr(iClust)
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
% Modified Dec 2004 JCJ -- Allows for multiple cluster inputs and shows which cluster is being operated on

redraw = false; rekey = false; undoable = false; % don't need to update

%%%% JCJ Dec 2004 - shows cluster being operated on %%%%
prompt={['Cluster(s) to Compare to cluster ' num2str(iClust) ':'],'bin size (msec):','Window width (msec):'};
def={'','1','500'};
dlgTitle='X Corr';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def);
if ~isempty(answer)
    %%%% JCJ Dec 2004 - processes multiple clusters in answer{1} vs. input cluster %%%% 
    iClust2=str2num(answer{1});
    for iC=1:length(iClust2)
        figure;
        MClustXcorr(iClust, iClust2(iC),  str2num(answer{2}),  str2num(answer{3}));
    end
end
