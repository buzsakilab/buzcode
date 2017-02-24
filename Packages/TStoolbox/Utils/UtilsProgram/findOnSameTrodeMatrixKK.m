function on_same_tet =   findOnSameTrodeMatrixKK(cell_names)
% on_same_tet =   findOnSameTet(cell_names)
%
% returns a matrix indicating the cell pairs belonging to the
% same trode. 
% INPUTS:
% cell_names = a cell array of strings containing the t-files names
%              (the TT  format is assumed)
% OUTPUT:
% on_same_tet = a matrix, elements being non-null for indices indicating
% cells on same trode

% fpbatta 2006
% status: developing



for i = 1:length(cell_names)
  idxs = findstr(cell_names{i}, '_');
  ci = cell_names{i};  
  c{i} = ci(1:idxs);
end

for i = 1:length(cell_names)
  for j = 1:length(cell_names)
    o(i,j) = strcmp(c{i}, c{j});
  end
end

on_same_tet = o;
  