function out = FindDuration(indx);

%FINDDURATION
%   Usage: out = FindDuration(indx)
%
% This is a short little function that finds a jump in the indx vector
% and then spits out the number of elements before that jump.  It does
% this for all jumps it finds.  By a jump, I mean a discontinutity in the
% index.  For example, if the index vector is [1 2 3 4 5 8 9 10], there is
% a jump between 5 and 8.

%Written by Dan Valente
%October 2007

if isempty(indx)
    out = [];
    return;
end

m = 1;
list = [];
indx(end+1)=-1;  %Makes sure last section of indx (past final juump) is properly considered.
for i = 1:length(indx)-1;
    if indx(i)+1 == indx(i+1)
        list = [list i];
    else
        list = [list i];
        out(m) = length(list);
        list = [];
        m = m+1;
    end
end




        
    