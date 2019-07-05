function bz_EditSessionMetadata(basepath)
% Sets up to allow user to edit session meta data in the specified session
% basepath.  Looks for a metadata file, if it doesn't find it, it creates a
% copy of the default template in the basepath.  Then opens the basepath
% .SessionMetadataText.m file in the editor for user input.  
% 
% INPUTS
%   basepath - computer path to the session folder of interest
%
% OUTPUTS
%   (Allows user to save a [basename].SessionMetadataText.m file in the basepath 
%
% Brendon Watson 2017

%% Initial variable parsing
if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end
basename = bz_BasenameFromBasepath(basepath);


%% Value setting - humans must do this.  Will be asked to edit a _NoteText.m file
notesname = [basename,'_SessionMetadataText.m'];
notesfullpath = fullfile(basepath,notesname);
if ~exist(notesfullpath,'file')
    w = which('bz_SessionMetadataTextTemplate.m');% copy an example header here to edit
    copyfile(w,notesname);
end

edit(notesname)

% prompt = 'Push any key in this window when done editing the NoteText file ';
% str = input(prompt,'s');
% run(notesname);%save _AnimalNotes.mat to disk
% 
% load(fullfile(basepath,[basename '.AnimalNotes.mat']))%load AnimalNotes
% AnimalMetadata = AnimalNotes;
% clear AnimalNotes