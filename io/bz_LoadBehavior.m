function [ behavior ] = bz_LoadBehavior( behaviorName,basePath )
%[ states ] = bz_LoadStates( statesName,baseName,datasetPath ) function for
%loading states.mat files.
%
%
%DLevenstein 2017
%%
if nargin < 2
    basePath = pwd;
end

baseName = bz_BasenameFromBasepath(basePath);
behaviorfile = fullfile(basePath,[baseName,'.',behaviorName,'.behavior.mat']);
evalin('caller',['load(''',behaviorfile,''')']);


end

