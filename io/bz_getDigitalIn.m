
function [digitalIn] = bz_getDigitalIn(ch,varargin)
% [pul, val, dur] = getPulses(d,varargin)
%
% Find digital In pulses
%
% INPUTS
% ch            Default all.
% <OPTIONALS>
% fs            Sampling frequency (in Hz), default 30000, or try to
%               recover for rhd 
% offset        Offset subtracted (in seconds), default 0.
% periodLag     How long a pulse has to be far from other pulses to be consider a different stimulation period (in seconds, default 5s)    
% filename      File to get pulses from. Default, digitalin.dat file with folder
%               name in current directory
%
%
% OUTPUTS
%               digitalIn - events struct with the following fields
% ints          C x 2  matrix with pulse times in seconds. First column of C 
%               are the beggining of the pulses, second column of C are the end of 
%               the pulses.
% dur           Duration of the pulses. Note that default fs is 30000.
% timestampsOn  Beggining of all ON pulses
% timestampsOff Beggining of all OFF pulses
% intsPeriods   Stimulation periods, as defined by perioLag
% 
% MV-BuzsakiLab 2019
% Based on Process_IntanDigitalChannels by P Petersen

% Parse options
if exist('ch') ~= 1
    ch = 'all';
end

p = inputParser;
addParameter(p,'fs',30000,@isnumeric)
addParameter(p,'offset',0,@isnumeric)
addParameter(p,'filename',[],@isstring)
addParameter(p,'periodLag',5,@isnumeric)

parse(p, varargin{:});
fs = p.Results.fs;
offset = p.Results.offset;
filename = p.Results.filename;
lag = p.Results.periodLag;

if ~isempty(dir('*.xml'))
    sess = bz_getSessionInfo(pwd,'noPrompts',true);
end
if ~isempty(dir('*DigitalIn.events.mat'))
    disp('Pulses already detected! Loading file.');
    file = dir('*DigitalIn.events.mat');
    load(file.name);
    return
end

if isempty(filename)
    filename=dir('digitalIn.dat');
    filename = filename.name;
else
    disp('No digitalIn file found...');
end

try [amplifier_channels, notes, aux_input_channels, spike_triggers,...
    board_dig_in_channels, supply_voltage_channels, frequency_parameters, board_adc_channels ] =...
    read_Intan_RHD2000_file_2;
    fs = frequency_parameters.board_dig_in_sample_rate;
catch
    disp('File ''info.rhd'' not found. (Type ''help <a href="matlab:help loadAnalog">loadAnalog</a>'' for details) ');
end

disp('Loading digital channels...');
m = memmapfile(filename,'Format','uint16','writable',false);
digital_word2 = double(m.Data);
clear m
Nchan = 16;
Nchan2 = 17;
for k = 1:Nchan
    tester(:,Nchan2-k) = (digital_word2 - 2^(Nchan-k))>=0;
    digital_word2 = digital_word2 - tester(:,Nchan2-k)*2^(Nchan-k);
    test = tester(:,Nchan2-k) == 1;
    test2 = diff(test);
    pulses{Nchan2-k} = find(test2 == 1);
    pulses2{Nchan2-k} = find(test2 == -1);
    data(k,:) = test;
end
digital_on = pulses;
digital_off = pulses2;
disp('Done!');

for ii = 1:size(digital_on,2)
    if ~isempty(digital_on{ii})
        % take timestamp in seconds
        digitalIn.timestampsOn{ii} = digital_on{ii}/fs;
        digitalIn.timestampsOff{ii} = digital_off{ii}/fs;
        
        % intervals
        d = zeros(2,max([size(digitalIn.timestampsOn{ii},1) size(digitalIn.timestampsOff{ii},1)]));
        d(1,1:size(digitalIn.timestampsOn{ii},1)) = digitalIn.timestampsOn{ii};
        d(2,1:size(digitalIn.timestampsOff{ii},1)) = digitalIn.timestampsOff{ii};
        if d(1,1) > d(2,1)
            d = flip(d,1);
        end
        if d(2,end) == 0; d(2,end) = nan; end
        digitalIn.ints{ii} = d;
        digitalIn.dur{ii} = digitalIn.ints{ii}(2,:) - digitalIn.ints{ii}(1,:); % durantion
        
        clear intsPeriods
        intsPeriods(1,1) = d(1,1); % find stimulation intervals
        intPeaks =find(diff(d(1,:))>lag);
        for jj = 1:length(intPeaks)
            intsPeriods(jj,2) = d(2,intPeaks(jj));
            intsPeriods(jj+1,1) = d(1,intPeaks(jj)+1);
        end
        intsPeriods(end,2) = d(2,end);  
        digitalIn.intsPeriods{ii} = intsPeriods;
    end
end

if exist('digitalIn')==1
    xt = linspace(0,size(data,2)/fs,size(data,2));
    data = flip(data);
    data = data(1:size(digitalIn.intsPeriods,2),:);

    h=figure;
    imagesc(xt,1:size(data,2),data);
    xlabel('s'); ylabel('Channels'); colormap gray 
    mkdir('Pulses');
    saveas(h,'pulses\digitalIn.png');

    try save([sess.FileName '.DigitalIn.events.mat'],'digitalIn');
    catch
        save('digitalIn.events.mat','digitalIn');
    end
else
    digitalIn = [];
end
 
end