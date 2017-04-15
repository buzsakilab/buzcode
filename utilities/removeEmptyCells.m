function [cell ind] = removeEmptyCells(cell)
ind = ~cellfun('isempty',cell);
cell = cell(~cellfun('isempty',cell));
end