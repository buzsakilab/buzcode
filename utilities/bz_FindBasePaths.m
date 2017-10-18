function [ basePaths,baseNames ] = bz_FindBasePaths(topPath)
%basePaths = bz_FindBasePaths(topPath) finds all buzcode basePaths in any
%directory within topPath and returns them.
%
%OUTPUT
%   basePaths   a cell array of folders that meet one the criteria for
%               buzcode basePaths
%
%Criteria
%   -has file baseName.xml
%   -has file baseName.sessionInfo.mat
%   ....where basePath is whateverPath/baseName/
%
%Please add more criteria to fit your needs! (Lines >31)
%DLevenstein 2017
%% Get all the directories
[basePaths, dirNames, fileNames] = dirwalk(topPath);

%% Loop through and check if theyre a buzcode basePath
for pp = length(basePaths):-1:1 %Count down to maintain indexing
    [ basepathYN ] = isBasePath(basePaths{pp},fileNames{pp});
    if ~basepathYN
        basePaths(pp) = [];
    end
end
baseNames = cellfun(@bz_BasenameFromBasepath,basePaths,'UniformOutput',false);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ basepathYN ] = isBasePath(pathName,thefiles)
%Given a candidate basePath (pathName) and the set of files in the path,
%this function returns yes if the folder fits criteria for buzcode basePath
thebasename = bz_BasenameFromBasepath(pathName);

xmlname = fullfile(pathName,[thebasename,'.xml']);
sessioninfoname = fullfile(pathName,[thebasename,'sessionInfo.mat']);

if exist(xmlname,'file') || exist(sessioninfoname,'file')
    basepathYN = true;
else
    basepathYN = false;
end

end