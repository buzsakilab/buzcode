function bz_PreprocessExtracellEphysSession(basepath)

if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end

%% Assuiming one already did bz_SetAnimalMetadata
bz_RunSessionMetadata(basepath);

%% making an xml if needed
basename = bz_BasenameFromBasepath(basepath);
txmlname = fullfile(basepath,[basename '.xml']);
d = dir(fullfile(basepath,txmlname));
load(fullfile(basepath,[basename '.SessionMetadata.mat']))
if isempty(d);%if no basename.xml already
    if isempty(SessionMetadata.ExtracellEphys.ParametersDivergentFromAnimalMetadata)%if not different from animal xml, copy the animal xml here
        ax = fullfile(SessionMetadata.AnimalMetadata.AnimalBasepath,[SessionMetadata.AnimalMetadata.AnimalName, '.xml']);
        copyfile(ax,txmlname);
    else    %if changes, make new xml for this session
        bz_MakeXML(basepath)
    end
end

%% 
disp('Concatenating .dat files')
bz_ConcatenateDats(basepath);
disp('Converting .dat to .lfp')
bz_LFPFromDat(basepath);