function bz_RunSessionMetadata(basepath)
% Wrapper to edit and then run SessionMetaData creation.
% This executes the metadata script created in bz_EditSessionMetadata.m 
% after it is created and the user hits a key to indicate as much.  
%
% See bz_EditSessionMetadata for much more info.  T
%
% Brendon Watson 2017

%% Initial variable parsing
if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end
basename = bz_BasenameFromBasepath(basepath);

%% 
bz_EditSessionMetadata(basepath)

%%
notesname = [basename,'_SessionMetadataText'];
notesfullpath = fullfile(basepath,notesname);
prompt = 'Push any key in this window when done editing the NoteText file ';
str = input(prompt,'s');
run(notesname);

