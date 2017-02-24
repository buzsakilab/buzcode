function [Y, I] = SortTable(X, L)



sx = X(:,L);

[sx, I] = sort(sx);

Y = X(I,:);

