function [sessionInfo] = bz_getSessionInfo(basePath)

if ~exist('basePath','var')
    basePath = pwd;
end

d = dir('*sessionInfo*');
if ~isempty(d) & length(d) == 1
   load(d.name); 
elseif ~isempty(d) & length(d) > 1
   warning('multiple sessionInfo files found, running LoadParameters instead..') 
   sessionInfo = LoadParameters(basePath);
else
   warning('could not find sessionInfo file, running LoadParameters instead..') 
   sessionInfo = LoadParameters(basePath);
end

%Here: check sessionInfo using bz_isSessionInfo

%Here: prompt user to add any missing sessionInfo fields and save

%Here: prompt user to save basePath/baseName.sessionInfo.mat if loaded from
%xml

end