function bz_ConcatenateDats(basepath,deleteoriginaldatsbool)
% bz_ConcatenateDats - Concatenate raw .dat files found in a session folder
% - for intan type recordings 
% 
% ALGORITHM OUTLINE: looks for .dat files in a folder (or in subfolders) to
% concatenate together.  The concatenation happens via system commands 
% ("cat" command for linux/mac, "copy" command if windows/pc).  Uses
% different assumptions to find and recognize relevant .dats depending on
% the acquisition system.  
% 
% REQUIREMENTS: Assumes you are in or pointed to a directory containing 
% subdirectories for various recording files from a single session. *It is 
% assumed that an earlier-acquired data file/folder will have a name that
% is sorted alphanumerically earlier.  Alphanumeric sorting order is
% assumed to be the recording temporal sequence.
% Works with acquisition systems: Intan  - 
%   1) intan: wherein subfolders are inside the session folder.  Each
%   subfolder contains simultaneously-recorded .dat files recorded for a
%   continuous period of time.  Start/stop recording commands each create a
%   new folder.  *It is assumed that the alphanumeric sorting of these 
%   folders corresponds with their sequence in acquisiton time.*  
%   These folders contain
%       - info.rhd files with metadata about the recording. 
%       - amplifier.dat - int16 file with usually neural data from the
%           headstage
%       - auxiliary.dat (optional) - uint16 file from auxiliary channels on
%           the headstage - often accelerometer
%       - analogin.dat (optional) - uint16 file from analogin channels on 
%           main board 
%       - digitalin.dat (optional) - uint16 file recording all 16 digital 
%           channels on main board 
%       - time.dat - int32 file giving recording sample indicies (e.g. 
%           0,1,2,3...) for each sample recorded in other channels
%       - supply.dat - uint16 file showing voltage supplied to (?preamp?)
%   
%
%  USAGE
%
%    bz_ConcatenateDats(basepath,deletedats)
%
%  INPUTS
%
%    basepath          computer path to session folder.  Defaults to
%                      current folder if no input given
%    deleteoriginaldatsbool  - boolean denoting whether to delete (1) or
%                              not delete (0) original .dats after
%                              concatenation.  Default = 0.
%
%  OUTPUT
%     Operates on files in specified folder.  No output variable
%
%  EXAMPLES
%      Can be called directly or via bz_PreprocessExtracellEphysSession.m
%
% Copyright (C) 2017 by Brendon Watson



%% Handling inputs
% basic session name and and path
if ~exist('basepath','var')
    basepath = cd;
end
basename = bz_BasenameFromBasepath(basepath);

if ~exist('deleteoriginaldatsbool','var')
    deleteoriginaldatsbool = 0;
end

%% assume xml is present in the basepath... comment this out... could copy an amplifier.xml from below
%if no xml, put one in the main path.
% if ~exist(fullfile(basepath,[basename,'.xml']),'file')
%     d = dir(basepath);
%     for a = 1:length(d)
%         if d(a).isdir 
%             if length(d(a).name) >= length(basename)
%                 if strcmp(d(a).name(1:length(basename)-2),basename(1:end-2))
%                     if exist(fullfile(basepath,d(a).name,'amplifier.xml'))
%                         disp(['copying xml from subfolder ' fullfile(basepath,d(a).name)])
%                         copyfile(fullfile(basepath,d(a).name,'amplifier.xml'),fullfile(basepath,[basename '.xml']))
%                         break
%                     end
%                 end
%             end
%         end
%     end
% end

%% If the dats are already merged quit
if exist(fullfile(basepath,[basename,'.dat']),'file')
    disp('.dat already exists in session directory, not merging subdats')
    return
end

%% Find all amplifier.dat paths and cat them to basename.dat in the session folder
otherdattypes = {'analogin';'digitalin';'auxiliary';'time';'supply'};
bad_otherdattypes = [];
for odidx = 1:length(otherdattypes)
    eval(['new' otherdattypes{odidx} 'path = fullfile(basepath,''' otherdattypes{odidx} '.dat'');'])
end

d = dir(basepath);
datpaths = {};
recordingbytes = [];
recordingnames = {};
rcount = 0;
for a = 1:length(d)
    if d(a).isdir 
%         if length(d(a).name) >= length(basename)
%             if strcmp(d(a).name(1:length(basename)-2),basename(1:end-2))
                if exist(fullfile(basepath,d(a).name,'amplifier.dat'))
                    % for unclear reasons the command below does not
                    % work with full paths as written... the line below
                    % that is the alternative but requires one to be in
                    % the right supradirectory
%                         datpaths{end+1} = fullfile(dirpath,d(a).name,'amplifier.dat');
                    rcount = rcount+1;
                    datpaths{rcount} = fullfile(d(a).name,'amplifier.dat');
                    t = dir(fullfile(basepath,d(a).name,'amplifier.dat'));
                    recordingbytes(rcount) = t.bytes;
                    recordingnames{rcount} = d(a).name;

                    for odidx = 1:length(otherdattypes)%loop through other .dat types found here
                        eval([otherdattypes{odidx} 'datpaths{rcount} = fullfile(basepath,recordingnames{rcount},''' otherdattypes{odidx} '.dat'');'])
                        eval(['d2 = dir(' otherdattypes{odidx} 'datpaths{rcount});'])
                        if isempty(d2)
                            bad_otherdattypes(odidx) = 1;
                        else
                            eval([otherdattypes{odidx} 'datsizes(rcount) = d2(1).bytes;'])
                        end
                    end

                end
%             end
%         end
    end
end
otherdattypes(find(bad_otherdattypes)) = [];%if there weren't analogin or digitalin in some recording
if isempty(datpaths)
    disp('No .dats found in subfolders.  Exiting bz_ConcatenateDats.')
    return
end

    
%% Concatenate
%     cs = strjoin(datpaths);
%     catstring = ['! cat ', cs, ' > ',fullfile(basepath,[basename,'.dat'])];
% 
%     % for some reason have to cd to supradirectory 
%     origdir = cd;
%     cd (basepath)
%     eval([catstring])
%     cd (origdir)
newdatpath = fullfile(basepath,[basename,'.dat']);
if isunix
    cs = strjoin(datpaths);
    catstring = ['! cat ', cs, ' > ',newdatpath];
elseif ispc%As of 4/9/2017 - never tested
    if length(datpaths)>1
        for didx = 1:length(datpaths)-1;
            datpathsplus{didx} = [datpaths{didx} '+'];
        end
    else
        datpathsplus = datpaths;
    end
    cs = strjoin(datpathsplus);
    catstring = ['! copy /b ', cs, ' ',newdatpath];
end

% action
eval(catstring)%execute concatention
    
%% Check that size of resultant .dat is equal to the sum of the components
%     t = dir(fullfile(basepath,[basename,'.dat']));
%     if t.bytes ~= sum(recordingbytes)
%         error('dat size not right')
%         return
%     end
% 
%     save(fullfile(basepath,[basename '_DatInfo.mat']),'recordingbytes','recordingnames')
t = dir(newdatpath);
if t.bytes ~= sum(recordingbytes)
    error('New .dat size not right.  Exiting')
    return
else
    disp(['Primary .dats concatenated and size checked'])
end

%% Also concatenate the other .dats
for odidx = 1:length(otherdattypes)
    eval(['tdatpaths = ' otherdattypes{odidx} 'datpaths;']);
    eval(['tnewdatpath = new' otherdattypes{odidx} 'path;']);
    if isunix
        cs = strjoin(tdatpaths);
        catstring = ['! cat ', cs, ' > ',tnewdatpath];
    elseif ispc%As of 4/9/2017 - never tested
        if length(tdatpaths)>1
            for didx = 1:length(tdatpaths)-1;
                datpathsplus{didx} = [tdatpaths{didx} '+'];
            end
        else
            datpathsplus = tdatpaths;
        end
        cs = strjoin(datpathsplus);
        catstring = ['! copy /b ', cs, ' ',tnewdatpath];
    end

    eval(catstring)%execute concatenation

    % Check that size of resultant .dat is equal to the sum of the components
    t = dir(tnewdatpath);
    eval(['recordingbytes = ' otherdattypes{odidx} 'datsizes;'])
    if t.bytes ~= sum(recordingbytes)
        error(['New ' otherdattypes{odidx} '.dat size not right.  Exiting after .dats converted.  Not deleting'])
        deleteoriginaldatsbool = 0;
    else
        disp([otherdattypes{odidx} ' concatenated and size checked'])
    end
end

%% Delete original dats, if that option is chosen
if deleteoriginaldatsbool
    %for other .dats
    for odidx = 1:length(otherdattypes)
        eval(['tdatpaths = ' otherdattypes{odidx} 'datpaths;']);
        for didx = 1:length(tdatpaths)
            if isunix 
                eval(['! rm ' tdatpaths{didx}])
            elseif ispc%As of 4/9/2017 - never tested
                eval(['! del ' tdatpaths{didx}])
            end
        end
    end
    %for main .dat
    for didx = 1:length(datpaths)
        if isunix 
            eval(['! rm ' datpaths{didx}])
        elseif ispc%As of 4/9/2017 - never tested
            eval(['! del ' datpaths{didx}])
        end
    end
end





















%%%% OLD BRENDON VERSION ASSUMING SESSIONMETADATA.MAT %%% 
% 
% % bz_ConcatenateDats - Concatenate raw .dat files found in a session folder
% % - for either intan type systems or amplipex
% % 
% % ALGORITHM OUTLINE: looks for .dat files in a folder (or in subfolders) to
% % concatenate together.  The concatenation happens via system commands 
% % ("cat" command for linux/mac, "copy" command if windows/pc).  Uses
% % different assumptions to find and recognize relevant .dats depending on
% % the acquisition system.  
% % 
% % REQUIREMENTS: Assumes you are in or pointed to a directory containing 
% % subdirectories for various recording files from a single session. *It is 
% % assumed that an earlier-acquired data file/folder will have a name that
% % is sorted alphanumerically earlier.  Alphanumeric sorting order is
% % assumed to be the recording temporal sequence.
% % Works with acquisition systems: Intan and Amplipex - 
% %   1) intan: wherein subfolders are inside the session folder.  Each
% %   subfolder contains simultaneously-recorded .dat files recorded for a
% %   continuous period of time.  Start/stop recording commands each create a
% %   new folder.  *It is assumed that the alphanumeric sorting of these 
% %   folders corresponds with their sequence in acquisiton time.*  
% %   These folders contain
% %       - info.rhd files with metadata about the recording. 
% %       - amplifier.dat - int16 file with usually neural data from the
% %           headstage
% %       - auxiliary.dat (optional) - uint16 file from auxiliary channels on
% %           the headstage - often accelerometer
% %       - analogin.dat (optional) - uint16 file from analogin channels on 
% %           main board 
% %       - digitalin.dat (optional) - uint16 file recording all 16 digital 
% %           channels on main board 
% %       - time.dat - int32 file giving recording sample indicies (e.g. 
% %           0,1,2,3...) for each sample recorded in other channels
% %       - supply.dat - uint16 file showing voltage supplied to (?preamp?)
% %   2) Amplipex: wherein there are a series of sequientially named .dat
% %   files, where the first might be called [basename]-01.dat, the next
% %   [basename]-02.dat etc.  *It is assumed that the alphanumeric sorting of
% %   these folders corresponds with their sequence in acquisiton time.*
% %   There are not assumed to be any sub-dat files, these .dats are assumed
% %   to carry all information.
% %
% %  USAGE
% %
% %    bz_ConcatenateDats(basepath,deletedats)
% %
% %  INPUTS
% %
% %    basepath          computer path to session folder.  Defaults to
% %                      current folder if no input given
% %    deleteoriginaldatsbool  - boolean denoting whether to delete (1) or
% %                              not delete (0) original .dats after
% %                              concatenation.  Default = 0.
% %
% %  OUTPUT
% %     Operates on files in specified folder.  No output variable
% %
% %  EXAMPLES
% %      Can be called directly or via bz_PreprocessExtracellEphysSession.m
% %
% % Copyright (C) 2017 by Brendon Watson
% % 
% % 
% % %% Input and directory handling 
% % if ~exist('basepath','var')
% %     basepath = cd;
% % elseif isempty(basepath)
% %     basepath = cd;
% % end
% % 
% % basename = bz_BasenameFromBasepath(basepath);
% % 
% % if ~exist('deletedats','var')
% %     deleteoriginaldatsbool = 0;
% % end
% %     
% % %%
% load(fullfile(basepath,[basename,'.SessionMetadata.mat']));
% newdatpath = fullfile(basepath,[basename,'.dat']);
% switch SessionMetadata.ExtracellEphys.RecordingSystem
%     case 'Intan'
%         otherdattypes = {'analogin';'digitalin';'auxiliary';'time';'supply'};
%         bad_otherdattypes = [];
%         for odidx = 1:length(otherdattypes)
%             eval(['new' otherdattypes{odidx} 'path = fullfile(basepath,''' otherdattypes{odidx} '.dat'');'])
%         end
%         
%         for fidx = 1:length(SessionMetadata.ExtracellEphys.Files.Names)% go through each folder in basepath to look for various .dats
%             datpaths{fidx} = fullfile(basepath,SessionMetadata.ExtracellEphys.Files.Names{fidx},'amplifier.dat');%int16
%             %datbytes already saved in SessionMetadata
%             
%             for odidx = 1:length(otherdattypes)%loop through other .dat types found here
%                 eval([otherdattypes{odidx} 'datpaths{fidx} = fullfile(basepath,SessionMetadata.ExtracellEphys.Files.Names{fidx},''' otherdattypes{odidx} '.dat'');'])
%                 eval(['d = dir(' otherdattypes{odidx} 'datpaths{fidx});'])
%                 if isempty(d)
%                     bad_otherdattypes(odidx) = 1;
%                 else
%                     eval([otherdattypes{odidx} 'datsizes(fidx) = d(1).bytes;'])
%                 end
%             end
% % 
% %             analogindatpaths{fidx} = fullfile(basepath,SessionMetadata.ExtracellEphys.Files.Names{fidx},'analogin.dat');
% %             d = dir(timedatpaths{fidx});
% %             analogindatsizes(fidx) = d(1).bytes;
% % 
% %             auxiliarydatpaths{fidx} = fullfile(basepath,SessionMetadata.ExtracellEphys.Files.Names{fidx},'auxiliary.dat');
% %             d = dir(timedatpaths{fidx});
% %             auxiliarydatsizes(fidx) = d(1).bytes;
% % 
% %             timedatpaths{fidx} = fullfile(basepath,SessionMetadata.ExtracellEphys.Files.Names{fidx},'time.dat');%int32
% %             d = dir(timedatpaths{fidx});
% %             timedatsizes(fidx) = d(1).bytes;
% % 
% %             supplydatpaths{fidx} = fullfile(basepath,SessionMetadata.ExtracellEphys.Files.Names{fidx},'supply.dat');%uint16
% %             d = dir(timedatpaths{fidx});
% %             supplydatsizes(fidx) = d(1).bytes;
%         end
%         otherdattypes(find(bad_otherdattypes)) = [];
%     case 'Amplipex'%As of 4/9/2017 - never tested
%         for fidx = 1:length(SessionMetadata.ExtracellEphys.Files.Names)
%             datpaths{fidx} = fullfile(basepath,[SessionMetadata.ExtracellEphys.Files.Names{fidx} '.dat']);
%         end        
% end
% 
% %% Concatenate the main data files
% if isunix
%     cs = strjoin(datpaths);
%     catstring = ['! cat ', cs, ' > ',newdatpath];
% elseif ispc%As of 4/9/2017 - never tested
%     if length(datpaths)>1
%         for didx = 1:length(datpaths)-1;
%             datpathsplus{didx} = [datpaths{didx} '+'];
%         end
%     else
%         datpathsplus = datpaths;
%     end
%     cs = strjoin(datpathsplus);
%     catstring = ['! copy /b ', cs, ' ',newdatpath];
% end
% 
% eval(catstring)%execute concatention
% 
% % Check that size of resultant .dat is equal to the sum of the components
% t = dir(newdatpath);
% recordingbytes = SessionMetadata.ExtracellEphys.Files.Bytes;
% if t.bytes ~= sum(recordingbytes)
%     error('New .dat size not right.  Exiting')
%     return
% else
%     disp(['Primary .dats concatenated successfully'])
% end
% 
% %% if intan, also concatenate the other .dats
% if strcmp(SessionMetadata.ExtracellEphys.RecordingSystem,'Intan')
%     for odidx = 1:length(otherdattypes)
%         eval(['tdatpaths = ' otherdattypes{odidx} 'datpaths;']);
%         eval(['tnewdatpath = new' otherdattypes{odidx} 'path;']);
%         if isunix
%             cs = strjoin(tdatpaths);
%             catstring = ['! cat ', cs, ' > ',tnewdatpath];
%         elseif ispc%As of 4/9/2017 - never tested
%             if length(tdatpaths)>1
%                 for didx = 1:length(tdatpaths)-1;
%                     datpathsplus{didx} = [tdatpaths{didx} '+'];
%                 end
%             else
%                 datpathsplus = tdatpaths;
%             end
%             cs = strjoin(datpathsplus);
%             catstring = ['! copy /b ', cs, ' ',tnewdatpath];
%         end
%         
%         eval(catstring)%execute concatenation
%         
%         % Check that size of resultant .dat is equal to the sum of the components
%         t = dir(tnewdatpath);
%         eval(['recordingbytes = ' otherdattypes{odidx} 'datsizes;'])
%         if t.bytes ~= sum(recordingbytes)
%             error(['New ' otherdattypes{odidx} '.dat size not right.  Exiting after .dats converted.  Not deleting'])
%             deleteoriginaldatsbool = 0;
%         else
%             disp([otherdattypes{odidx} ' concatenated successfully'])
%         end
%     end
% end
% 
% 
% %% Delete original dats, if that option is chosen
% if deleteoriginaldatsbool
%     %for other .dats
%     for odidx = 1:length(otherdattypes)
%         eval(['tdatpaths = ' otherdattypes{odidx} 'datpaths;']);
%         for didx = 1:length(tdatpaths)
%             if isunix 
%                 eval(['! rm ' tdatpaths{didx}])
%             elseif ispc%As of 4/9/2017 - never tested
%                 eval(['! del ' tdatpaths{didx}])
%             end
%         end
%     end
%     %for main .dat
%     for didx = 1:length(datpaths)
%         if isunix 
%             eval(['! rm ' datpaths{didx}])
%         elseif ispc%As of 4/9/2017 - never tested
%             eval(['! del ' datpaths{didx}])
%         end
%     end
% end

