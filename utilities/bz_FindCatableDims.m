function [ catabledim ] = bz_FindCatableDims( tryarrays )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%INPUT
%   tryarrays   cell array of arrays you'd like to try to concatenate
%
%OUTPUT
%   catabledims dimensions along which you can concatenate
%
%DLevenstein 2018
%%
%Remove any empty arrays
tryarrays(cellfun(@(X) isempty(X),tryarrays)) = [];

arraysizes = cellfun(@size,tryarrays,'uniformoutput',false);
numdims = cellfun(@length,arraysizes);

%If all the arrays don't have the same number of dimensions, you're out of luck
if range(numdims)~=0
    disp('Your arrays don''t have the same number of dimensions. You''re SOL')
end

arraysizes = cat(1,arraysizes{:});

samesizes = range(arraysizes,1)==0;

if all(samesizes)
    catabledim = 1:numdims;
elseif sum(~samesizes)==1
    catabledim = find(~samesizes);
else
    disp('Two or more of your dimensions are different size. You''re SOL')
end


end

