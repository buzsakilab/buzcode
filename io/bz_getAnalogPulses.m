
function [pulses] = bz_getAnalogPulses(varargin)
% [pul, val, dur] = bz_getAnalogPulses(varargin)
%
% Find square pulses. If not argument it provided it tries to find pulses
% in intan analog in file.
%
% <OPTIONALS>
% analogCh      List of analog channels with pulses to be detected (it support Intan Buzsaki Edition).
% data          R x C matrix with analog data. C is data, R should be
%               greater than 1.
% fs            Sampling frequency (in Hz), default 30000.
% offset        Offset subtracted (in seconds), default 0.
% periodLag     How long a pulse has to be far from other pulses to be consider a different stimulation period (in seconds, default 5s)    
% filename      File to get pulses from. Default, data file with folder
%               name in current directory
% manualThr     Check manually threslhold amplitude (default, false)
%
%
% OUTPUTS
%               pulses - events struct with the following fields
% timestamps    C x 2  matrix with pulse times in seconds. First column of C 
%               are the beggining of the pulses, second column of C are the end of 
%               the pulses. 
% amplitude     values of the pulses with respect baleline (normalized as 0).
% duration      Duration of the pulses. Note that default fs is 30000.
% eventID       Numeric ID for classifying various event types (C X 1)
% eventIDlabels label for classifying various event types defined in eventID (cell array C X 1)  
% intsPeriods   Stimulation periods, as defined by perioLag
%
% Manu-BuzsakiLab 2018

% Parse options
p = inputParser;
addParameter(p,'analogCh',[],@isnumeric)
addParameter(p,'data',[],@isnumeric)
addParameter(p,'fs',30000,@isnumeric)
addParameter(p,'offset',0,@isnumeric)
addParameter(p,'filename',[],@isstring)
addParameter(p,'periodLag',5,@isnumeric)
addParameter(p,'manualThr',false,@islogical) %% to do!!!

parse(p, varargin{:});
fs = p.Results.fs;
offset = p.Results.offset;
filename = p.Results.filename;
lag = p.Results.periodLag;
manualThr = p.Results.manualThr;
d = p.Results.data;
analogCh = p.Results.analogCh;

sess = bz_getSessionInfo(pwd,'noPrompts',true);
if exist([sess.FileName '.pulses.events.mat'],'file') 
    disp('Pulses already detected! Loading file.');
    load([sess.FileName '.pulses.events.mat']);
    return
end

if isempty(d) && isempty(filename)                                         % is d is not a signal, and not filename specified
    disp('No filename... looking for pulses file...');
    f = [];
    if exist('analogin.dat','file') == 2
        f=dir('analogin.dat');
    end
    if isempty(f) || f.bytes == 0                                              % if analogin is empty or doesn't exist
        disp('analogin.dat file is empty or does not exist, was the recording made in Intan Buzsaki edition?');
        f = dir('*amplifier*.dat');                                        % is exist amplifier
        if isempty(f)
            filename = split(pwd,'\'); filename = filename{end};
            filename = strcat(filename,'.dat');
        else
            filename = f.name;
        end
        nCh = sess.nChannels;
        disp('Loading file...');
        tic
        d = bz_LoadBinary(filename, 'frequency', fs, 'nChannels', nCh,'channels',analogCh+1);
        toc
    else
        disp('Loading analogin.dat...');
        filename = 'analogin.dat';  
        try [amplifier_channels, notes, aux_input_channels, spike_triggers,...         
            board_dig_in_channels, supply_voltage_channels, frequency_parameters, board_adc_channels ] =...
            read_Intan_RHD2000_file_2;
        catch
            disp('File ''info.rhd'' not found. (Type ''help <a href="matlab:help loadAnalog">loadAnalog</a>'' for details) ');
        end
        num_channels = length(board_adc_channels); % ADC input info from header file
        fileinfo = dir(filename);
        disp('Loading file...');
        num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes
        fid = fopen('analogin.dat', 'r');
        v = fread(fid, [num_channels, num_samples], 'uint16');
        fclose(fid);
        d = v(analogCh,:);
        d = (d - 32768) * 0.0003125; % convert to volts
        fs = frequency_parameters.amplifier_sample_rate;
    end
end

if size(d,1) > size(d,2)
    d = d';
end

for jj = 1 : size(d,1)
    xt = linspace(1,length(d)/fs,length(d));
    fprintf(' ** Channel %3.i of %3.i... \n',jj, size(d,1));
    if any(d<0) % if signal go negative, rectify
        d = d - min(d);
    end
    
    if ~manualThr
        thr = 25*median(abs(d(jj,:))/0.6745); % computing threshold
        if thr == 0 || ~any(d>thr)
            disp('Trying 5*std threshold...');
            d = d - mean(d);
            thr = 5 * std(double(d(jj,:)));
        end
    else
        h = figure;
        plot(xt(1:100:end), d(1:100:end));
        xlabel('s'); ylabel('amp');
        title('Select threshold with the mouse and press right click...');
        [~,thr] = ginput(1);
        hold on
        plot([xt(1) xt(end)],[thr thr],'-r');
        pause(1);
        close(h);
    end
    
    dBin = (d(jj,:)>thr); % binarize signal
    locsA = find(diff(dBin)==1)/fs; % start of pulses
    locsB = find(diff(dBin)==-1)/fs; % end of pulses

    pul{jj} = locsA(1:min([length(locsA) length(locsB)]));
    for ii = 1 : size(pul{jj},2) % pair begining and end of the pulse
        try pul{jj}(2,ii) =  locsB(find(locsB - pul{jj}(1,ii) ==...
            min(locsB(locsB > pul{jj}(1,ii)) - pul{jj}(1,ii))));
        catch
            keyboard;
        end
    end

    baseline_d = int32(median(d(jj,:)));
    val{jj}=[];
    for ii = 1 : size(pul{jj},2) % value of the pulse respect basaline
        val{jj}(ii) = median(int32(d(jj,int32(pul{jj}(1,ii) * fs : pul{jj}(2,ii) * fs)))) - baseline_d;
    end
    
    pul{jj} = pul{jj} - offset;
    % discard pulses < 2 * median(abs(x)/0.6745) as noise or pulses in negatives times
    idx = find((val{jj} < thr*0.4) | pul{jj}(1,:)<0);
    val{jj}(idx) = [];
    pul{jj}(:,idx) = [];
    
    if ~isempty(pul{jj})
        dur{jj} = pul{jj}(2,:) - pul{jj}(1,:); % durantion
        
        stimPer{jj}(1,1) = pul{jj}(1,1); % find stimulation intervals
        intPeaks =find(diff(pul{1}(1,:))>lag);
        for ii = 1:length(intPeaks)
            stimPer{jj}(ii,2) = pul{jj}(2,intPeaks(ii));
            stimPer{jj}(ii+1,1) = pul{jj}(1,intPeaks(ii)+1);
        end
        stimPer{jj}(end,2) = pul{jj}(2,end);
    else
        dur{jj} = [];
        stimPer{jj} = [];
    end
    
    eventID{jj} = ones(size(dur{jj})) * jj;
    
    h=figure;
    subplot(1,size(d,1),jj);
    hold on
    plot(xt(1:100:end), d(1:100:end));
    plot(xt([1 end]), [thr thr],'r','LineWidth',2);
    xlim([0 xt(end)]);
    ax = axis;
    if ~isempty(locsA)
        plot(locsA, ax(4),'o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','none','MarkerSize',3);
    end
    xlabel('s'); ylabel('Amplitude (au)'); 
end
mkdir('Pulses');
saveas(gca,'pulses\pulThr.png');

if ~isempty(locsA) % if no pulses, not save anything...
    pulses.timestamps = cell2mat(pul)';
    pulses.amplitude = cell2mat(val)';
    pulses.duration = cell2mat(dur)';
    pulses.intsPeriods = cell2mat(stimPer);
    pulses.eventID = cell2mat(eventID)';
    disp('Saving locally...');
    save([sess.FileName '.pulses.events.mat'],'pulses');
else
    pulses = [];
end
 
end