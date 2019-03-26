function [ sessionInfo,success ] = bz_tagChannel( basePath,channums,tag,varargin )
%[ sessionInfo,success ] = bz_tagChannel( basePath,channums,tag ) tag the channels in channums with
%the tag (string).
%Will be saved in the sessionInfo as sessionInfo.channelTags.tag
%
%(options)
%   'noPrompts'
%   'overwrite'     (default: false), adds channel to tag instead of overwriting
%DLevenstein 2019
%%
p = inputParser;
addParameter(p,'noPrompts',false,@islogical);
addParameter(p,'overwrite',false,@islogical);

parse(p,varargin{:})
noPrompts = p.Results.noPrompts;
overwrite = p.Results.overwrite;

%%
%Load the sessionInfo
[sessionInfo] = bz_getSessionInfo(basePath,'noPrompts',noPrompts);

%Tag the channels
if ~isfield(sessionInfo,'channelTags') || ~isfield(sessionInfo.channelTags,tag) || overwrite
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

