function bz_ConcatenateDats(basepath,deletedats)
%assumes you are in or pointed to a directory containing subdirectories for
% various recording files from a single session

%% Input and directory handling 
if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end

basename = bz_BasenameFromBasepath(basepath);

if ~exist('deletedats','var')
    deletedats = 0;
end
    
%%
load(fullfile(basepath,[basename,'.SessionMetadata.mat']));
newdatpath = fullfile(basepath,[basename,'.dat']);
switch SessionMetadata.ExtracellEphys.RecordingSystem
    case 'Intan'
        otherdattypes = {'analogin';'digitalin';'auxiliary';'time';'supply'};
        bad_otherdattypes = [];
        for odidx = 1:length(otherdattypes)
            eval(['new' otherdattypes{odidx} 'path = fullfile(basepath,''' otherdattypes{odidx} '.dat'');'])
        end
        
        for fidx = 1:length(SessionMetadata.ExtracellEphys.Files.Names)% go through each folder in basepath to look for various .dats
            datpaths{fidx} = fullfile(basepath,SessionMetadata.ExtracellEphys.Files.Names{fidx},'amplifier.dat');%int16
            %datbytes already saved in SessionMetadata
            
            for odidx = 1:length(otherdattypes)%loop through other .dat types found here
                eval([otherdattypes{odidx} 'datpaths{fidx} = fullfile(basepath,SessionMetadata.ExtracellEphys.Files.Names{fidx},''' otherdattypes{odidx} '.dat'');'])
                eval(['d = dir(' otherdattypes{odidx} 'datpaths{fidx});'])
                if isempty(d)
                    bad_otherdattypes(odidx) = 1;
                else
                    eval([otherdattypes{odidx} 'datsizes(fidx) = d(1).bytes;'])
                end
            end
% 
%             analogindatpaths{fidx} = fullfile(basepath,SessionMetadata.ExtracellEphys.Files.Names{fidx},'analogin.dat');
%             d = dir(timedatpaths{fidx});
%             analogindatsizes(fidx) = d(1).bytes;
% 
%             auxiliarydatpaths{fidx} = fullfile(basepath,SessionMetadata.ExtracellEphys.Files.Names{fidx},'auxiliary.dat');
%             d = dir(timedatpaths{fidx});
%             auxiliarydatsizes(fidx) = d(1).bytes;
% 
%             timedatpaths{fidx} = fullfile(basepath,SessionMetadata.ExtracellEphys.Files.Names{fidx},'time.dat');%int32
%             d = dir(timedatpaths{fidx});
%             timedatsizes(fidx) = d(1).bytes;
% 
%             supplydatpaths{fidx} = fullfile(basepath,SessionMetadata.ExtracellEphys.Files.Names{fidx},'supply.dat');%uint16
%             d = dir(timedatpaths{fidx});
%             supplydatsizes(fidx) = d(1).bytes;
        end
        otherdattypes(find(bad_otherdattypes)) = [];
    case 'Amplipex'%As of 4/9/2017 - never tested
        for fidx = 1:length(SessionMetadata.ExtracellEphys.Files.Names)
            datpaths{fidx} = fullfile(basepath,[SessionMetadata.ExtracellEphys.Files.Names{fidx} '.dat']);
        end        
end

%% Concatenate the main data files
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

eval(catstring)%execute concatention

% Check that size of resultant .dat is equal to the sum of the components
t = dir(newdatpath);
recordingbytes = SessionMetadata.ExtracellEphys.Files.Bytes;
if t.bytes ~= sum(recordingbytes)
    error('New .dat size not right.  Exiting')
    return
else
    disp(['Primary .dats concatenated successfully'])
end

%% if intan, also concatenate the other .dats
if strcmp(SessionMetadata.ExtracellEphys.RecordingSystem,'Intan')
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
            deletedats = 0;
        else
            disp([otherdattypes{odidx} ' concatenated successfully'])
        end
    end
end


%% Delete original dats, if that option is chosen
if deletedats
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

