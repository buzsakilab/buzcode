function [inds_all, inds_cell] = intvl_to_inds(Nx2_intvls, N_ONLY)
% inds = intvl_to_inds(Nx2_intvls, N_ONLY)
%
% Created by EWS, Oct 2012

if (nargin < 2)
    N_ONLY = 0;
end
inds_cell = cell(size(Nx2_intvls,1),1);
for i=1:size(Nx2_intvls,1)
    inds_cell{i} = (Nx2_intvls(i,1):Nx2_intvls(i,2))';
end
inds_all = unique(cell2mat(inds_cell));
if N_ONLY
    inds_all = length(unique(cell2mat(inds_cell)));
    inds_cell = cellfun(@length, inds_cell);
end
