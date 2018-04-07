function bz_RunAnimalMetadata(basepath)

%% Initial variable parsing
if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end
basename = bz_BasenameFromBasepath(basepath);

%% 
bz_EditAnimalMetadata(basepath)

%%
notesname = [basename,'_AnimalMetadataText'];
notesfullpath = fullfile(basepath,notesname);
prompt = 'Push any key in this window when done editing the NoteText file ';
str = input(prompt,'s');
run(notesfullpath);
