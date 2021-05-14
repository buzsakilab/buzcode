function bz_PreprocessExtracellEphysAnimal(basepath)

if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end

%% Assuiming one already did bz_SetAnimalMetadata
bz_RunAnimalMetadata(basepath);

%% making an xml if needed
bz_MakeXML(basepath)
