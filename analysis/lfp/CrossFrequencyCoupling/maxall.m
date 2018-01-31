function [maxval, max_subs] = maxall(X)
% [maxval, max_subs] = maxall(X)
%
% Return the maximum value and the subscripts of the maximum value in X.
%
% Created by EW Schomburg, April 2013

[maxval, tmpind] = max(X(:));
if (length(size(X)) == 2) && (nnz(size(X) == 1) == 1)
    max_subs = tmpind;
else
    [tmp_subs{1:length(size(X))}] = ind2sub(size(X), tmpind);
    max_subs = cell2vec(tmp_subs);
    max_subs = max_subs(:)';
end
