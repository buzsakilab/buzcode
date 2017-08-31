function [sessionInfo] = bz_getSessionInfo(basePath)

if ~exist('basePath','var')
    basePath = pwd;
end

d = dir('*sessionInfo*');
if ~isempty(d) & length(d) == 1
   load(d.name); 
else
   warning('could not find sessionInfo file, running LoadParameters instead..') 
   sessionInfo = LoadParameters(basePath);
end