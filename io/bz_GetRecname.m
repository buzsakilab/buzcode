function [ recName ] = bz_GetRecname( basePath )
%Returns the recording name of a basePath formatted like:
% whateverFolder/recName/
%%
[datasetfolder,recName] = fileparts(basePath);


end

