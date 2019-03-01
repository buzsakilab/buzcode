function [  ] = bz_Counter( current,total,name )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if current==1
    fprintf('%s (of %d): %d',name,total,current)
elseif current==total
    for j=0:log10(current-1)
      fprintf('\b'); % delete previous counter display
    end
    fprintf('\b DONE!\n')
else
    for j=0:log10(current-1)
      fprintf('\b'); % delete previous counter display
    end
    fprintf('%d',current)
end


end

