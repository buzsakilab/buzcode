function [  ] = bz_Counter(current,total,name)
%bz_Counter(current,total,name) counts up loop indices, given current index
%and the total number of indices. Displays a name of what's being counted.

%Add.... capability or fix for parforloop...

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

