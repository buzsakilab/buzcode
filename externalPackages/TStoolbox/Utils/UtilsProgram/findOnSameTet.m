function on_same_tet =   findOnSameTet(cell_names)
% on_same_tet =   findOnSameTet(cell_names)

% returns a vector of element indicating the cell pairs belonging to the
% same tetrode. 
% INPUTS:
% cell_names = a cell array of strings containing the t-files names
%              (the BBClust filename format is assumed)
% OUTPUT:
% on_same_tet = a vector suitable as an argument for Stephen's Multi_partial

% batta 2001
% status: developing



for i = 1:length(cell_names)
  c{i} = strtok(cell_names{i}, '_');
end

for i = 1:length(cell_names)
  for j = 1:length(cell_names)
    o(i,j) = strcmp(c{i}, c{j});
  end
end

on_same_tet = find(o);
  