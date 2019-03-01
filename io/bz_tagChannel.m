function [ sessionInfo,success ] = bz_tagChannel( basePath,channums,tag )
%bz_tagChannel( basePath,channums,tag ) tag the channels in channums with
%the tag (string).
%Will be saved in the sessionInfo as sessionInfo.channelTags.tag
%
%DLevenstein 2019
%%

%Load the sessionInfo
[sessionInfo] = bz_getSessionInfo(basePath);

%Tag the channels
if ~isfield(sessionInfo,'channelTags') || ~isfield(sessionInfo.channelTags,tag)
    sessionInfo.channelTags.(tag) = channums;
else %If the tag alredy exists, add these channels to that tag
    sessionInfo.channelTags.(tag) = unique([sessionInfo.channelTags.(tag) channums]);
end

%Save the sessionInfo
baseName = bz_BasenameFromBasepath(basePath);
filename = fullfile(basePath,[baseName,'.sessionInfo.mat']);
try
    save(filename,'sessionInfo');
    success=true;
catch
    success=false;
    display('Failed to save sessionInfo with tag');
end

end

