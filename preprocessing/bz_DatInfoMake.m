function DatInfo = bz_DatInfoMake(basepath)
% Store basic info about the original dat files: names and bytes
% Brendon Watson 2016-8

%% Input and directory handling 
if ~exist('basepath','var')
    basepath = cd;
elseif isempty(basepath)
    basepath = cd;
end

basename = bz_BasenameFromBasepath(basepath);


%% First find out what system was used to record this... based on files present
% check if amplipex (.metas present) or intan (rhd.infos present)... or neither
recsys = '';

%Amplipex case
dmeta = dir(fullfile(basepath,'*.meta'));
if ~isempty(dmeta)
    recsys = 'Amplipex';
end

%Intan case
drhd = dir(basepath);
for didx = 1:length(drhd)
   if drhd(didx).isdir
       t = dir(fullfile(basepath,drhd(didx).name,'info.rhd'));
       if ~isempty(t)%only use folders with info.rhd inside
            recsys = 'Intan';
            break
       end
   end
end

%% Now gather info
switch recsys
    case 'Amplipex'%if from an amplipex
        d = dmeta;
        DatInfo.RecordingSystem = 'Amplipex';
        for idx = 1:length(d);
            DatInfo.Files.Names{idx} = d(idx).name(1:end-5);
            tdat = strcat(d(idx).name(1:end-5),'.dat');
            t = dir(fullfile(basepath,tdat));
            if ~isempty(t) %if original .dat still around
                DatInfo.Files.Bytes(idx) = t.bytes;
            else%if no .dat around (ie anymore), get bytes from .meta 
                disp([tdat ' not found']);
                DatInfo.Files.Bytes{idx} = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'filebytes'));
            end
            NChannels(idx) = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'NumChans'));
            SamplingRate(idx) = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'SamplingRate'));
            Amplification(idx) = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'Gain'));
        end
        if ~isempty(NChannels) %if some files were found
            %Channels
            if sum(abs(diff(NChannels)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different numbers of channels - per .meta files.  Will use count from first recording.')
            end
            DatInfo.Parameters.NumberOfChannels = NChannels(1);

            %Sampling Rate
            if sum(abs(diff(SamplingRate)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different SamplingRates - per .meta files.  Will use number from first recording.')
            end
            DatInfo.Parameters.SampleRate = SamplingRate(1);

            %Amplification
            if sum(abs(diff(Amplification)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different Amplifications - per .meta files.  Will use number from first recording.')
            end
            DatInfo.Parameters.Amplification = Amplification(1);

            DatInfo.Parameters.VoltsPerUnit = VoltsPerUnit_Amplirec(basename,basepath);%calculate from .metas/.inis
            DatInfo.Parameters.BitsPerSample = 16;%amplirec default
            DatInfo.Parameters.VoltageRange = 10;%not used except to make xml
        end
    case 'Intan' %if not amplipex, look for intan
        d = drhd;
        NAmpChannels = [];
        NAuxChannels = [];
        NDigitalChannels = [];
        SamplingRate = [];

        for didx = 1:length(d)
            if d(didx).isdir
                t = dir(fullfile(basepath,d(didx).name,'info.rhd'));
                if ~isempty(t)%only use folders with info.rhd inside

                   DatInfo.RecordingSystem = 'Intan';
                   rhd = fullfile(basepath,d(didx).name,'info.rhd');
    %                ampdats{end+1} = fullfile(basepath,d(didx).name,'amplifier.dat');
    %                recordingdiridxs = didx;
    
                   if ~isfield(DatInfo,'Files')% if first runthrough
                       DatInfo.Files.Names = [];
                       DatInfo.Files.Bytes = [];
                   end
                   
                   DatInfo.Files.Names{end+1} = d(didx).name;
                   
                   % read the info.rhd file for each individual recorded
                   % file/folder
                   try
                      [amplifier_channels, notes, aux_input_channels, spike_triggers,...         
                          board_dig_in_channels, supply_voltage_channels, frequency_parameters ] =...
                          read_Intan_RHD2000_file(fullfile(basepath,d(didx).name),'info.rhd');%function from Intan
                      NAmpChannels(end+1) = length(amplifier_channels);
                      NAuxChannels(end+1) = length(aux_input_channels);
                      NDigitalChannels(end+1) = length(board_dig_in_channels);
                      SamplingRate(end+1) = frequency_parameters.amplifier_sample_rate;

                   catch%if read error, set values based on prior
                       if ~isempty(NAmpChannels)
%                            if ~isnan(NAmpChannels(1));
                              NAmpChannels(end+1) = NAmpChannels(1);
%                            else
%                               NAmpChannels(end+1) = SessionMetadata.AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal;
%                            end
                       else
                           NAmpChannels(end+1) = nan;
                       end
                       if ~isempty(NAuxChannels)
                           NAuxChannels(end+1) = NAuxChannels(1);
                       else
                           NAuxChannels(end+1) = nan;
                       end
                       if ~isempty(NDigitalChannels)
                           NDigitalChannels(end+1) = NDigitalChannels(1);
                       else
                           NDigitalChannels(end+1) = nan;
                       end
                       if ~isempty(SamplingRate)
                           SamplingRate(end+1) = SamplingRate(1);
                       else
                           SamplingRate(end+1) = nan;
                       end
                   end

                   t = dir(fullfile(basepath,d(didx).name,'amplifier.dat'));
                   if ~isempty(t) %if original .dat still around
                       DatInfo.Files.Bytes(end+1) = t.bytes;
                   else%if no .dat around
                       disp([DatInfo.Files.Names{end} ' not found']);
                       timedat = dir(fullfile(basepath,d(didx).name,'time.dat'));
                       if ~isnan(NAmpChannels(end));
                          DatInfo.Files.Bytes(end+1) = timedat(1).bytes/2*NAmpChannels(end);
                       else
                          DatInfo.Files.Bytes(end+1) = timedat(1).bytes/2*SessionMetadata.AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal;
                       end
                   end

               end
            end
        end
        %handle errors created by improper read of rhd.info files
        if isnan(sum(abs(diff(NAmpChannels))))
            warning('At least one read error on rhd.info files.  Using values from other files to populate.')
            ok = find(~isnan(NAmpChannels),1,'first');
            if isempty(ok)
                warning('No good reads on NumberofAmpChannels. Assuming based on probe maps')
                NAmpChannels(isnan(NAmpChannels)) = SessionMetadata.AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal
            else
                NAmpChannels(isnan(NAmpChannels)) = NAmpChannels(ok);
            end
        end
        if isnan(sum(abs(diff(NAuxChannels))))
            warning('At least one read error on rhd.info files.  Using values from other files to populate.')
            ok = find(~isnan(NAuxChannels),1,'first');
            if isempty(ok)
                error('No good reads on NAuxChannels, assuming 0')
                NAuxChannels(isnan(NAuxChannels)) = 0;
            else
                NAuxChannels(isnan(NAuxChannels)) = NAuxChannels(ok);
            end
        end
        if isnan(sum(abs(diff(NDigitalChannels))))
            warning('At least one read error on rhd.info files.  Using values from other files to populate.')
            ok = find(~isnan(NDigitalChannels),1,'first');
            if isempty(ok)
                error('No good reads on NDigitalChannels, assuming 0')
                NDigitalChannels(isnan(NDigitalChannels)) = 0;
            else
                NDigitalChannels(isnan(NDigitalChannels)) = NDigitalChannels(ok);
            end
        end
        if isnan(sum(abs(diff(SamplingRate))))
            warning('At least one read error on rhd.info files.  Using values from other files to populate.')
            ok = find(~isnan(SamplingRate),1,'first');
            if isempty(ok)
                error('No good reads on SamplingRate, taking from AnimalMetadata')
                SamplingRate(isnan(SamplingRate)) = SessionMetadata.AnimalMetadata.ExtracellEphys.Parameters.SampleRate;
            else
                SamplingRate(isnan(SamplingRate)) = SamplingRate(ok);
            end
        end

        %handle errors where different parameters read for different files (rather than uniform)            
        if sum(abs(diff(NAmpChannels)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different numbers of channels - per info.rhd files.  Will use count from first recording.')
        end
        DatInfo.Parameters.NumberOfChannels = NAmpChannels(1);

        if sum(abs(diff(NAuxChannels)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different numbers auxiliary of channels - per info.rhd files.  Will use count from first recording.')
        end
        DatInfo.NumAuxiliaryChannels = NAuxChannels(1);

        if sum(abs(diff(NDigitalChannels)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different numbers digital of channels - per info.rhd files.  Will use count from first recording.')
        end
        DatInfo.NumDigitalChannels = NAuxChannels(1);

        %sampling rate
        if sum(abs(diff(SamplingRate)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different SamplingRates - per info.rhd files.  Will use number from first recording.')
        end
        DatInfo.Parameters.SampleRate = SamplingRate(1);

        % Not gather-able from intan info.rhd files
        DatInfo.Parameters.Amplification = 1;%digitized on chip, let's say 1
        DatInfo.Parameters.VoltsPerUnit = 0.0000002;%intan default                
        DatInfo.Parameters.BitsPerSample = 16;%intan default
        DatInfo.Parameters.VoltageRange = 10;%not used except to make xml

        %Now that everything is in order, calculate seconds from bytes
        DatInfo.Files.Seconds = DatInfo.Files.Bytes/DatInfo.Parameters.NumberOfChannels/2/DatInfo.Parameters.SampleRate;        
    case '' %if no known recording system used
        warning('No recording system found, using defaults.  May manually set values in guts of this .m file')
        DatInfo.RecordingSystem = 'Unknown';
        DatInfo.Files.Names = {};
        DatInfo.Files.Bytes = [];
        DatInfo.Files.Seconds = [];%can manually write in here if desired

        DatInfo.Parameters.NumberOfChannels = SessionMetadata.AnimalMetadata.ExtracellEphys.Channels.NumChannelsTotal;

        DatInfo.Parameters.SampleRate = 20000;%Lab default
        DatInfo.Parameters.Amplification = 1;%Intan digitized on chip, let's say 1
        DatInfo.Parameters.VoltsPerUnit = 0.0000002;%Intan default                
        DatInfo.Parameters.BitsPerSample = 16;%Intan default
        DatInfo.Parameters.VoltageRange = 10;%not used except to make xml
end
    
save(fullfile(basepath,[basename '_DatInfo.mat']),'DatInfo')
