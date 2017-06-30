
function [vals inds equal_max_val_inds] = randmax(X)
% The function is the same as the regular max function, but if there are ties in maximum value, 
%   this function returns a randomly selected best index rather than always choosing 
%   the index with the smallest value (as would happen if normal max is used).  This function
%   can return a variable equal_max_values which is 1 if there are two or more entries 
%   that are equal to the maximum value, and 0 otherwise.  

    

% if the X is a row vector take its transpose to make it a column vector
if size(X, 1) == 1
    X = X';
end

equal_max_values = 0;


if size(size(X), 2) > 2
    error('max rand only capable right now of dealing with vectors and 2D matrices')
end



[vals inds] = max(X);
[sorted_vals sorted_inds] = sort(X, 'descend');

% if there are two equal max values, choose one randomly
equal_max_val_inds = sorted_vals(1, :) == sorted_vals(2, :);   % if this is a vector and there are ties, this will equal 1, and the below code still works (b/c iTies is 1)

if sum(equal_max_val_inds) > 0
    
    for iTies = find(equal_max_val_inds)   % go through each column of the matrix that has tied values and...
         curr_max_inds = sorted_inds(find(sorted_vals(:, iTies) == sorted_vals(1, iTies)), iTies);   % find all indecies that have the maximum value
         curr_max_inds = curr_max_inds(randperm(length(curr_max_inds)));   % randomly select an index from all the ones that have the maximum value       
         inds(iTies) = curr_max_inds(1);        
    end
    
end







