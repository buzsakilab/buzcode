function bz_MakeXML(basepath)
% Starts from a basepath to make a .xml file.  Assumes presence of
% .SessionMetadata.mat and .AnimalMetadata.mat
%
% Brendon Watson 2017

if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end
basename = bz_BasenameFromBasepath(basepath);

%%
if exist(fullfile(basepath,[basename '.SessionMetadata.mat']))
    load(fullfile(basepath,[basename '.SessionMetadata.mat']))
    params = SessionMetadata.ExtracellEphys.Parameters;
    pfiles = SessionMetadata.AnimalMetadata.ExtracellEphys.Probes.ProbeLayoutFilenames;
    plugord = SessionMetadata.ExtracellEphys.Probes.PluggingOrder;
elseif exist(fullfile(basepath,[basename '.AnimalMetadata.mat']))
    load(fullfile(basepath,[basename '.AnimalMetadata.mat']))
    pfiles = AnimalMetadata.ExtracellEphys.Probes.ProbeLayoutFilenames;
    plugord = AnimalMetadata.ExtracellEphys.Probes.PluggingOrder;
    params = AnimalMetadata.ExtracellEphys.Parameters;
end
    
%% Make an initial .xml file
bz_MakeXMLFromProbeMaps(basepath,basename,pfiles,plugord,params);
