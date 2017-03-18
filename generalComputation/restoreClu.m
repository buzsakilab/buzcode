function [] = restoreClu(basepath,basename,shank)
try
if ~isempty(dir('*kwik'))
    if nargin < 1
    basepath = pwd;
    name = dir('*alg*');
    s = split(name.name,'.');
    basename = s{1};
    shank = str2num(s{end});
    end
    
    try
        tkwik = fullfile(basepath,num2str(shank),[basename '_sh' num2str(shank) '.kwik']);


        clu_orig = h5read(tkwik,['/channel_groups/' num2str(shank) '/spikes/clusters/original']);
    catch
        tkwik = fullfile(basepath,[basename '_sh' num2str(shank) '.kwik']);


        clu_orig = h5read(tkwik,['/channel_groups/' num2str(shank) '/spikes/clusters/original']);
    end
h5write(tkwik,['/channel_groups/' num2str(shank) '/spikes/clusters/main'],uint32(clu_orig));
disp(['restored ' tkwik])
end
catch
    
   disp('couldnt restore, corrupt file?') 
end
end
