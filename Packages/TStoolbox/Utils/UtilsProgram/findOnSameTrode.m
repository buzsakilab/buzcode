function on_same_tet =   findOnSameTrode(cell_names)
% on_same_tet =   findOnSameTet(cell_names)

% returns a vector of element indicating the cell pairs belonging to the
% same Warp trode. 
% INPUTS:
% cell_names = a cell array of strings containing the t-files names
%              (the BBClust filename format is assumed)
% OUTPUT:
% on_same_tet = a vector suitable as an argument for Stephen's Multi_partial

% batta 2001
% status: developing



for i = 1:length(cell_names)
  idxs = findstr(cell_names{i}, '_');
  ci = cell_names{i};  
  c{i} = ci(1:idxs(end)-1);
end

for i = 1:length(cell_names)
  for j = 1:length(cell_names)
    o(i,j) = strcmp(c{i}, c{j});
  end
end

on_same_tet = find(o);
  