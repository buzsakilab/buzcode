function DatsMetadata = bz_DatFileMetadata(basepath)
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
for didx = length(drhd):-1:1
   ok = 0;
   if drhd(didx).isdir
       t = dir(fullfile(basepath,drhd(didx).name,'info.rhd'));
       if ~isempty(t)%only use folders with info.rhd inside
            recsys = 'Intan';
            ok = 1;
       end
   end
   if ~ok
        drhd(didx) = [];
   end
end

%% Now gather info
switch recsys
    case 'Amplipex'%if from an amplipex
        d = dmeta;
        DatsMetadata.RecordingSystem = 'Amplipex';
        for idx = 1:length(d);
            DatsMetadata.Recordings.Names{idx} = d(idx).name(1:end-5);
            tdat = strcat(d(idx).name(1:end-5),'.dat');
            t = dir(fullfile(basepath,tdat));
            if ~isempty(t) %if original .dat still around
                DatsMetadata.Recordings.Bytes(idx) = t.bytes;
            else%if no .dat around (ie anymore), get bytes from .meta 
                disp([tdat ' not found']);
                DatsMetadata.Recordings.Bytes{idx} = str2num(bz_ReadAmplipexMetafileAspects(fullfile(basepath,d(idx).name),'filebytes'));
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
            DatsMetadata.Parameters.NumberOfChannels = NChannels(1);

            %Sampling Rate
            if sum(abs(diff(SamplingRate)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different SamplingRates - per .meta files.  Will use number from first recording.')
            end
            DatsMetadata.Parameters.SampleRate = SamplingRate(1);

            %Amplification
            if sum(abs(diff(Amplification)))%if different numbers of channels differed across recordings in session, error out
                warning('Rercordings contain different Amplifications - per .meta files.  Will use number from first recording.')
            end
            DatsMetadata.Parameters.Amplification = Amplification(1);

            DatsMetadata.Parameters.VoltsPerUnit = VoltsPerUnit_Amplirec(basename,basepath);%calculate from .metas/.inis
            DatsMetadata.Parameters.BitsPerSample = 16;%amplirec default
            DatsMetadata.Parameters.VoltageRange = 10;%not used except to make xml
        end
    case 'Intan' %if not amplipex, look for intan
        d = drhd;
        NAmpChannels = [];
        NAuxChannels = [];
        NDigitalChannels = [];
        SamplingRate = [];
        thisidx = 0; 
        
        for didx = 1:length(d)
            if d(didx).isdir
                t = dir(fullfile(basepath,d(didx).name,'info.rhd'));
                if ~isempty(t)%only use folders with info.rhd inside

                   DatsMetadata.RecordingSystem = 'Intan';
                   rhd = fullfile(basepath,d(didx).name,'info.rhd');
    %                ampdats{end+1} = fullfile(basepath,d(didx).name,'amplifier.dat');
    %                recordingdiridxs = didx;
    
                   thisidx = thisidx + 1;
                  
                   DatsMetadata.Recordings.Names{thisidx} = d(didx).name;
                   
                   % read the info.rhd file for each individual recorded
                   % file/folder
                   try
                      [amplifier_channels, notes, aux_input_channels, spike_triggers,...         
                          board_dig_in_channels, supply_voltage_channels, frequency_parameters ] =...
                          read_Intan_RHD2000_file(fullfile(basepath,d(didx).name),'info.rhd');%function from Intan
                      DatsMetadata.Recordings.IntanRHDInfoRaw{thisidx} = v2struct(amplifier_channels, notes, aux_input_channels, spike_triggers,...         
                          board_dig_in_channels, supply_voltage_channels, frequency_parameters);
                      NAmpChannels(thisidx) = length(amplifier_channels);
                      NAuxChannels(thisidx) = length(aux_input_channels);
                      NDigitalChannels(thisidx) = length(board_dig_in_channels);
                      SamplingRate(thisidx) = frequency_parameters.amplifier_sample_rate;

                   catch%if read error, set values based on prior
                       if ~isempty(NAmpChannels)
%                            if ~isnan(NAmpChannels(1));
                          NAmpChannels(thisidx) = NAmpChannels(1);
                       else
                          sess = bz_getSessionInfo(basepath,'noPrompts',true);
                          NAmpChannels(thisidx) = sess.nChannels;
                       end
                       if ~isempty(NAuxChannels)
                           NAuxChannels(thisidx) = NAuxChannels(1);
                       else
                           NAuxChannels(thisidx) = nan;
                       end
                       if ~isempty(NDigitalChannels)
                           NDigitalChannels(thisidx) = NDigitalChannels(1);
                       else
                           NDigitalChannels(thisidx) = nan;
                       end
                       if ~isempty(SamplingRate)
                           SamplingRate(thisidx) = SamplingRate(1);
                       else
                           SamplingRate(thisidx) = nan;
                       end
                   end

                   t = dir(fullfile(basepath,d(didx).name,'amplifier.dat'));
                   t2 = dir(fullfile(basepath,d(didx).name,'amplifier_analogin_auxiliary_int16.dat'));
                   if ~isempty(t2) %if original .dat still around.  Handling Luke's RHD software to put all 16bit into same file
                       DatsMetadata.Recordings.Bytes(thisidx) = t2.bytes;                       
                   elseif ~isempty(t) %if original .dat still around
                       DatsMetadata.Recordings.Bytes(thisidx) = t.bytes;
                   else%if no .dat around
                       disp([DatsMetadata.Recordings.Names{end} ' not found']);
                       timedat = dir(fullfile(basepath,d(didx).name,'time.dat'));
                       if ~isnan(NAmpChannels(end));
                          DatsMetadata.Recordings.Bytes(thisidx) = timedat(1).bytes/2*NAmpChannels(end);
                       else
                           sess = bz_getSessionInfo(basepath);
                           DatsMetadata.Recordings.Bytes(thisidx) = timedat(1).bytes/2*sess.nChannels;
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
                sess = bz_getSessionInfo(basepath);
                NAmpChannels(isnan(NAmpChannels)) = sess.nChannels;
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
                sess = bz_getSessionInfo(basepath);
                SamplingRate(isnan(SamplingRate)) = sess.rates.wideband;
            else
                SamplingRate(isnan(SamplingRate)) = SamplingRate(ok);
            end
        end

        %handle errors where different parameters read for different files (rather than uniform)            
        if sum(abs(diff(NAmpChannels)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different numbers of channels - per info.rhd files.  Will use count from first recording.')
        end
        DatsMetadata.Parameters.NumberOfChannels = NAmpChannels(1);

        if sum(abs(diff(NAuxChannels)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different numbers auxiliary of channels - per info.rhd files.  Will use count from first recording.')
        end
        DatsMetadata.NumAuxiliaryChannels = NAuxChannels(1);

        if sum(abs(diff(NDigitalChannels)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different numbers digital of channels - per info.rhd files.  Will use count from first recording.')
        end
        DatsMetadata.NumDigitalChannels = NAuxChannels(1);

        %sampling rate
        if sum(abs(diff(SamplingRate)))%if different numbers of channels differed across recordings in session, error out
            warning('Rercordings contain different SamplingRates - per info.rhd files.  Will use number from first recording.')
        end
        DatsMetadata.Parameters.SampleRate = SamplingRate(1);

        % Not gather-able from intan info.rhd files
        DatsMetadata.Parameters.Amplification = 1;%digitized on chip, let's say 1
        DatsMetadata.Parameters.VoltsPerUnit = 0.0000002;%intan default                
        DatsMetadata.Parameters.BitsPerSample = 16;%intan default
        DatsMetadata.Parameters.VoltageRange = 10;%not used except to make xml

        %Now that everything is in order, calculate seconds from bytes
        DatsMetadata.Recordings.Seconds = DatsMetadata.Recordings.Bytes/DatsMetadata.Parameters.NumberOfChannels/2/DatsMetadata.Parameters.SampleRate;        
    case '' %if no known recording system used
        warning('No recording system found, using defaults.  May manually set values in guts of this .m file')
        DatsMetadata.RecordingSystem = 'Unknown';
        DatsMetadata.Recordings.Names = {};
        DatsMetadata.Recordings.Bytes = [];
        DatsMetadata.Recordings.Seconds = [];%can manually write in here if desired

        try 
            sess = bz_getSessionInfo(basename);
        catch
            disp('No metadata or session info detected, quitting bz_DatFileMetadata')
            return
        end
        
        DatsMetadata.Parameters.NumberOfChannels = sess.nChannels;
        DatsMetadata.Parameters.SampleRate = sess.rates.wideband;%Lab default
        DatsMetadata.Parameters.Amplification = 1;%Intan digitized on chip, let's say 1
        DatsMetadata.Parameters.VoltsPerUnit = 0.0000002;%Intan default                
        DatsMetadata.Parameters.BitsPerSample = sess.nBits;%Intan default
        DatsMetadata.Parameters.VoltageRange = 10;%not used except to make xml
end
    
save(fullfile(basepath,[basename '_DatsMetadata.mat']),'DatsMetadata')
