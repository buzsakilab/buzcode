function bz_PreprocessExtracellEphysSession(basepath)

if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end

%% Assuming one already did bz_SetAnimalMetadata
bz_RunSessionMetadata(basepath);

%% making an xml if needed
basename = bz_BasenameFromBasepath(basepath);
txmlname = fullfile(basepath,[basename '.xml']);
d = dir(txmlname);
load(fullfile(basepath,[basename '.SessionMetadata.mat']))
if isempty(d);%if no basename.xml already
    if isempty(SessionMetadata.ExtracellEphys.ParametersDivergentFromAnimalMetadata)%if not different from animal xml, copy the animal xml here
        ax = fullfile(SessionMetadata.AnimalMetadata.AnimalBasepath,[SessionMetadata.AnimalMetadata.AnimalName, '.xml']);
        if exist(ax,'file')
            copyfile(ax,txmlname);
        end
    end
end
d = dir(txmlname);
if isempty(d);%if no basename.xml already
    bz_MakeXML(basepath)
end

%% 
disp('Concatenating .dat files')
bz_ConcatenateDats(basepath);
disp('Converting .dat to .lfp')
bz_LFPFromDat(basepath);