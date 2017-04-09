function bz_PreprocessExtraEphysSession(basepath,basename)

if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end

if ~exist('basename','var')
    [~,basename] = fileparts(basepath);
elseif isempty(basename)
    [~,basename] = fileparts(basepath);
end


bz_SetSessionMetadata(basepath,basename);
bz_ConcatenateDats(basepath,basename);
bz_LFPFromDat(basepath,basename);