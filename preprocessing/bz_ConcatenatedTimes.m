function  bz_ConcatenatedTimes(basepath,sortFiles)
%   bz_ConcatenatedTimes(basepath)
%   This function generates the events.mat file that saves all the merge
%   information. Useful if you dind't do this when you merged the dats.
%
%  USAGE
%
%    bz_ConcatenateDats(basepath,sortFiles)
%
%  INPUTS
%
%    basepath          computer path to session folder.  Defaults to
%                      current folder if no input given
%    sortFiles               - boolean denoting whether to sort files according 
%                              to time of recording (1) or
%                              not (0) and thus sort them alphabetically 
%                              Default = 1.
%
%  OUTPUT
%     Operates on files in specified folder.  No output variable
%
%   Antonio FR 10/2018


%% Handling inputs
if ~exist('basepath','var')
    basepath = cd;
end
basename = bz_BasenameFromBasepath(basepath);

if ~exist('sortFiles','var')
    sortFiles = 1;
end

%% If the dats are not merged quit
if ~exist(fullfile(basepath,[basename,'.dat']),'file')
    disp('no merged dat in directory, quitting')
    return
end

if exist(fullfile(basepath,[basename,'.MergePoints.events.mat']),'file')
    disp('MergedPoints already in directory')
    return
end

%% Find all .dat paths in subfolders 

otherdattypes = {'analogin';'digitalin';'auxiliary';'time';'supply'};
bad_otherdattypes = [];
for odidx = 1:length(otherdattypes)
    %eval(['new' otherdattypes{odidx} 'path = fullfile(basepath,''' otherdattypes{odidx} '.dat'');'])
    newpaths.(otherdattypes{odidx}) = fullfile(basepath,[otherdattypes{odidx}, '.dat']);
end

d = dir(basepath);
datpaths = {};
datsizes.amplifier = [];
recordingnames = {};
rcount = 0; %Count of good subfolders
for a = 1:length(d)
    %look in each subfolder
    if d(a).isdir 
        %Check for amplifier.dat or subfolderbaseName.dat 
        if exist(fullfile(basepath,d(a).name,[d(a).name,'.dat']),'file')
            ampfile = fullfile(basepath,d(a).name,[d(a).name,'.dat']);
        else
            ampfile = fullfile(basepath,d(a).name,'amplifier.dat');
        end
        
        if exist(ampfile,'file')
            rcount = rcount+1;
            datpaths.amplifier{rcount} = ampfile;
            t = dir(ampfile);
            datsizes.amplifier(rcount) = t.bytes;
            recordingnames{rcount} = d(a).name;

            for odidx = 1:length(otherdattypes)%loop through other .dat types found here
               % eval([otherdattypes{odidx} 'datpaths{rcount} = fullfile(basepath,recordingnames{rcount},''' otherdattypes{odidx} '.dat'');'])
                datpaths.(otherdattypes{odidx}){rcount} = fullfile(basepath,recordingnames{rcount},[otherdattypes{odidx} '.dat']);
                %eval(['d2 = dir(' otherdattypes{odidx} 'datpaths.amplifier{rcount});'])
                d2 = dir(datpaths.(otherdattypes{odidx}){rcount});
                if isempty(d2)
                    bad_otherdattypes(odidx) = 1;
                else
                    %eval([otherdattypes{odidx} 'datsizes(rcount) = d2(1).bytes;'])
                    datsizes.(otherdattypes{odidx})(rcount) = d2(1).bytes;
                end
            end
        end
    end
end
otherdattypes(find(bad_otherdattypes)) = [];%if there weren't analogin or digitalin in some recording
if isempty(datpaths.amplifier)
    disp('No .dats found in subfolders.  Exiting bz_ConcatenateDats.')
    return
end

%% Get the XML
try 
    %Look for xml/sessionInfo in topfolder
    sessionInfo = bz_getSessionInfo(basepath,'noPrompts',true);
catch
    %If none exists, look for xml in any of the subpaths
    disp('No .xml or .sessionInfo in top folder, trying subfolders')
    for ff = 1:length(recordingnames)
        try
            sessionInfo = LoadParameters(fullfile(basepath,recordingnames{ff}));
            xmlfilename = fullfile(sessionInfo.session.path,[sessionInfo.session.name,'.xml']);
            [SUCCESS,MESSAGE,MESSAGEID] = copyfile(...
                xmlfilename,...
                fullfile(basepath,[basename,'.xml']),'f');
            display(['Copied xml from ',recordingnames{ff}])
            break
        catch
        end
    end
end

%% Sort files according to time of recording

if sortFiles

    try
        names2sort = cellfun(@(X) str2num(X(end-5:end)),recordingnames,'UniformOutput',false);
        names2sort = cell2mat(names2sort);
        %disp('Assuming the last 6 digits reflect recording time.')
        %disp('Don''t like it? Write in some new options for sorting.')
    catch
        disp('Last 6 digits not numeric... sorting alphabetically')
    end

    [~,I] = sort(names2sort);
    recordingnames = recordingnames(I);
    datpaths.amplifier = datpaths.amplifier(I);
    datsizes.amplifier = datsizes.amplifier(I);
    for odidx = 1:length(otherdattypes)
        datpaths.(otherdattypes{odidx}) = datpaths.(otherdattypes{odidx})(I);
        datsizes.(otherdattypes{odidx}) = datsizes.(otherdattypes{odidx})(I);
    end
else
    disp('not implemented for non temporally sorted concatenation')
end

%% Check that size of resultant .dat is equal to the sum of the components
newdatpath = fullfile(basepath,[basename,'.dat']);
t = dir(newdatpath);
if t.bytes ~= sum(datsizes.amplifier)
    error('New .dat size not right.  Exiting')
    return
else
    sizecheck.amplifier = true;
    disp('Primary .dats concatenated and size checked')
end

%% Get time points from the time.dat
%Use the timestamps from time.dat to get the sort order
%Number of samples in time.dat. First timepoint, last timepoint
%Convert from number of samples to recording time of start/ends
for ff = 1:length(datpaths.time)
    
	f = fopen(datpaths.time{ff},'r'); 
    % Determine total number of samples in file
    fileStart = ftell(f);
    
    %Read the first time point
    firsttimepoint = fread(f,1,'int32');
    status = fseek(f,-4,'eof'); %int32 = 4 bytes
    lasttimepoint = fread(f,1,'int32');
    fileStop = ftell(f);
    
    firstlasttimepoints(ff,:) = [firsttimepoint lasttimepoint];
    numsamples(ff) = fileStop./4;
    if ff==1
        transitiontimes_samp = firstlasttimepoints(ff,:);
    else
        transitiontimes_samp(ff,:) = firstlasttimepoints(ff,:)+transitiontimes_samp(ff-1,2)+1;
    end
end

disp(['Calculating merge times based on wideband samplingRate of ',num2str(sessionInfo.rates.wideband),'Hz.'])
transitiontimes_sec = transitiontimes_samp./sessionInfo.rates.wideband; %convert to seconds

%% Make the events.mat file that saves all the merge information
% filesT.int == MergePoint.timestamps

eventsfilename = fullfile(basepath,[basename,'.MergePoints.events.mat']);

MergePoints.timestamps = transitiontimes_sec;
MergePoints.timestamps_samples = transitiontimes_samp;
MergePoints.firstlasttimpoints_samples = firstlasttimepoints;
MergePoints.foldernames = recordingnames;
MergePoints.filesmerged = datpaths;
MergePoints.filesizes = datsizes;
MergePoints.sizecheck = sizecheck;
MergePoints.detectorinfo.detectorname = 'bz_ConcatenateDats';
MergePoints.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');


save(eventsfilename,'MergePoints');


end

