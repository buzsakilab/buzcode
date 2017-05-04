function bz_WriteEventFile(filepath,filename,eventTimes,eventLabels)

% bz_WriteEventFile(filepath,filename,SR, eventTimes)

% DESCRIPTION  This fx takes event times and generates a neuroscope .evt file

% INPUTS
% filepath:    Path where want to save .evt file
% filename:    Name of .evt file; if [], output .evt will have same name as .xml file in filepath
% eventTimes:  Matrix of event times in seconds, in format time by event x type
% eventLabels: Cell array of column header labels (strings) for each event type

% OUTPUT       .evt file in designated filepath

% E.G.         Save output of bz_FindRipples as a .evt file for viewing in neuroscope
%              bz_WriteEventFile(filepath,[],[ripples.times(:,1) ripples.peaks ripples.times(:,2)],{'tstart'; 'tpeak'; 'tstop'})

%% Event file name
% Ensure don't write over extant .evt file
evtFiles = dir([filepath '*.evt']);
if isempty(evtFiles)
    fileN = 1;
else
    %set file index to next available value\
    pat = '.R[0-9].';
    fileN = 0;
    for ii = 1:length(evtFiles)
        token  = regexp(evtFiles(ii).name,pat);
        val    = str2double(evtFiles(ii).name(token+2:token+4));
        fileN  = max([fileN val]);
    end
    fileN = fileN + 1;
end


if isempty(filename)
    datafile = dir([filepath '*.xml']);
    filename = datafile.name(1:end-4);
else
    error('Filename not provided and no .xml file in filepath. Must provide filename.')
end

%% Convert detections to milliseconds (neuroscope default)
eventTimes = eventTimes.*(1000);

%% Write event file
fid = fopen(sprintf('%s%s%s.R%02d.evt',filepath,filesep,filename,fileN),'w');
fprintf(1,'Writing event file ...\n');

for eventtime = 1:size(eventTimes,1)
    for eventtype = 1:size(eventTimes,2)
        fprintf(fid,['%9.1f\' eventLabels{eventtype} '\n'],eventTimes(eventtime,eventtype));
    end
end

fclose(fid);
