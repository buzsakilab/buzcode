function [ behavior ] = bz_ConcatenateBehavior( behaviorName, basePath, varargin)
%[ behavior ] = bz_ConcatenateBehavior( behaviorName, basePath)
%Concatenates buzcode behavior.mat files from subfolders, using the 
%MergePoints from bz_ConcatenateDats. Automatically adjusts timestamps from
%each of the subfiles to reflect the timestamps in the merged .lfp file.
%Saves a new file: basePath/baseName.behaviorName.behavior.mat
%
%INPUTS
%   behaviorName    assumes that one or more subfolders have a file called
%                   subfolder.behaviorName.behavior.mat, following buzcode
%                   behavior.mat guidelines.
%                   Automatically adjusts timestamps stored in
%                   behaviorName.timestamps and subfields of
%                   behaviorName.ints.
%   basePath        upper level basePath into which you have already run
%                   bz_ConcatenateDats. 
%                   Requires basePath/baseName.MergePoints.events.mat
%   (options)
%   'timefields'    any fields that contain timestamps that you'd like to
%                   adjust to match the concatenation MergePoints
%                   (note: this feature is not yet implemented, let me know
%                   if you need it and we'll add it to the function...)
%
%note: automatically adjusts behaviorName.timestamps and any fields in
%behaviorName.ints to account for merge points. Any additional fields can 
%be indicated with optional 'timefields' input
%DLevenstein 2019
%% input Parsing

p = inputParser;
addParameter(p,'timefields',{});
addParameter(p,'remerge',false);

parse(p,varargin{:})

timefields = p.Results.timefields;
REMERGE = p.Results.remerge;

%% DEV
% basePath = '/Users/dlevenstein/Desktop/180530_KO_EM1M3';
% behaviorName = 'pupildiameter';

%%
baseName = bz_BasenameFromBasepath(basePath);
savename = fullfile(basePath,[baseName,'.',behaviorName,'.behavior.mat']);

if exist(savename,'file') && ~REMERGE
    RESULT = questdlg([savename,' already exists, overwrite?  ',...
        '(Use ''remerge'',true to avoid this message in the future)'],...
        'Yes','No');
    switch RESULT
        case {'No','Cancel'}
            return
    end
end
    
%%
%Load the mergepoints metadata file
MergePoints = bz_LoadEvents(basePath,'MergePoints');
if isempty(MergePoints)
    error(['No MergePoints.events.mat file found in ',basePath])
end



%%

subfolders = cellfun(@(X) fullfile(basePath,X),MergePoints.foldernames,...
    'UniformOutput',false);
behaviorfilenames = cellfun(@(X,Y) fullfile(X,[Y,'.',behaviorName,'.behavior.mat']),...
    subfolders,MergePoints.foldernames,'UniformOutput',false);

%Check that all the subfolders, and which behavior files, exist
existsubfolders = cellfun(@(X) logical(exist(X,'dir')),subfolders);
existbehfiles = cellfun(@(X) logical(exist(X,'file')),behaviorfilenames);
existbehfiles = find(existbehfiles);


numsubfolders = length(subfolders);

if any(~existbehfiles)
    display(['Can''t find ',behaviorName,'.behavior.mat in ', MergePoints.foldernames(~existbehfiles),...
        '.  Merging from ',MergePoints.foldernames(existbehfiles)])
end

%%
offsets = MergePoints.timestamps(:,1);
subbeh = struct([]);
for ff = 1:length(existbehfiles)
	thisidx = existbehfiles(ff);
	newbeh = bz_LoadBehavior(subfolders{thisidx},behaviorName);

    % check that the fields match, if not add the new field 
    %(add notification here if not alinged)
    [subbeh,newbeh] = bz_Matchfields(subbeh,newbeh,'add');
	subbeh(ff) = newbeh;
    
	% Required fields: timestamps, data. Should prompt user for other
	% timestamp fields.
	subbeh(ff).timestamps = subbeh(ff).timestamps+offsets(thisidx);
    
    if isfield(subbeh(ff),'ints')
        intfields = fieldnames(subbeh(ff).ints);
        for gg = 1:length(intfields)    
            subbeh(ff).ints.(intfields{gg}) = subbeh(ff).ints.(intfields{gg})+offsets(thisidx);
        end
    end
    
end


%Collapse and concatenate the fields from the structures
behavior = bz_CollapseStruct(subbeh,'match','justcat',true);

tol = 1e-4; %tolerance for difference in sample rates
if range(behavior.samplingRate)>tol
    warning('Sampling Rates are not the same...... taking the average')
end
behavior.samplingRate = mean(behavior.samplingRate);

behavior.MergeInfo.MergedFiles = behaviorfilenames(existbehfiles);
behavior.MergeInfo.MergedDate = datetime('today');


%% Save the new behavior file
eval([behaviorName,'=behavior;'])
save(savename,behaviorName)

end

