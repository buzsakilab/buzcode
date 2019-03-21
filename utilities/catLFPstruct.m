function [lfpCat] = catLFPstruct(lfp)
% [lfpCat] = catLFPstruct(lfp)
%   
% Concatenate inervals from an buzcode lfp structure 

lfpCat = lfp(1);
lfpCat.data = []; lfpCat.timestamps = []; lfpCat.duration = []; lfpCat.interval = [];
for i = 1:size(lfp,2)
    lfpCat.data = cat(1,lfpCat.data,lfp(i).data);
    lfpCat.timestamps = cat(1,lfpCat.timestamps,lfp(i).timestamps);
    lfpCat.duration = cat(1,lfpCat.duration,lfp(i).duration);
    lfpCat.interval = cat(1,lfpCat.interval,lfp(i).interval);
end

end


