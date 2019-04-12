function CutDatFragment(fname,timeperiod,RenameOriginalsAsOrig)
% Creating a new dat file from a fragment of an original dat file
% INPUTs
% fname: filepath/name of original dat
% timeperiod: start and stop of file fragment to keep, in seconds
% RenameOriginalsAsOrig: [optional] - specifies naming convention for new vs
% old files.  
%    If RenameOriginalsAsOrig = 1 then amplifier.dat becomes
%      amplifier.dat_orig and the new file is called amplifier.dat.  Same with
%      supply.dat becoming supply.dat etc.
%    If RenameOriginalsAsOrig = 0 then amplifier.dat remains and the new
%      file is caled amplifier_fragment.dat.  Same for other files like
%      supply, time etc
%    Default = 0
%
% Original by Antonio Ruiz
% Modified by Peter Petersen, Brendon Watson
%
% function should be run from the directory of the dat files
% It assumes an info.rhd in the directory
%
% fname = 'Peter_MS13_171130_143907';
% timeperiod = [0,44*60+9];

%% Input handling
if ~exist('RenameOriginalsAsOrig','var')
    RenameOriginalsAsOrig = 0;
end

if strcmp(fname(end-3:end),'.dat')
    fname = fname(1:end-4);
end


%% Gather recording meta info
[amplifier_channels, notes, aux_input_channels, spike_triggers,...         
board_dig_in_channels, supply_voltage_channels, frequency_parameters ] = read_Intan_RHD2000_file(pwd,'info.rhd');

% amplifier_channels = 1:1:64;
% aux_input_channels = 3;
% supply_voltage_channels = 1;
% frequency_parameters.amplifier_sample_rate = 30000;

%% amplifier.dat
disp('Writing amplifier file')
NumCh = length(amplifier_channels);
SampRate = frequency_parameters.amplifier_sample_rate;
if RenameOriginalsAsOrig
    inname = [fname '.dat_orig'];
    outname = [fname '.dat'];
    movefile(outname,inname)
else
    inname = [fname '.dat'];
    outname = [fname '_fragment.dat'];
end
a = memmapfile(inname,'Format','int16');
aa = a.data;
fid = fopen(outname,'W');
for i = timeperiod(1)+1:timeperiod(2)
     fwrite(fid,aa((i-1)*NumCh*SampRate+1:i*NumCh*SampRate),'int16');
end
fclose(fid);
clear aa a

%% auxiliary.dat
disp('Writing aux file')
if RenameOriginalsAsOrig
    inname = 'auxiliary.dat_orig';
    outname = 'auxiliary.dat';
    movefile(outname,inname)
else
    inname = 'auxiliary.dat';
    outname = 'auxiliary_fragment.dat';
end
NumCh = length(aux_input_channels);
m = memmapfile(inname,'Format','uint16');
h1 = fopen(outname,'W');
for i = timeperiod(1)+1:timeperiod(2)
    fwrite(h1,m.Data((i-1)*NumCh*SampRate+1:i*NumCh*SampRate),'uint16');
end
fclose(h1);

%% analogin.dat
d = dir('analogin.dat');
if ~isempty(d)
    disp('Writing analogin file')
    if RenameOriginalsAsOrig
        inname = 'analogin.dat_orig';
        outname = 'analogin.dat';
        movefile(outname,inname)
    else
        inname = 'analogin.dat';
        outname = 'analogin_fragment.dat';
    end
    NumCh = length(board_adc_channels);
    m = memmapfile(inname,'Format','uint16');
    h2 = fopen(outname,'W');
    for i = timeperiod(1)+1:timeperiod(2)
        fwrite(h2,m.Data((i-1)*NumCh*SampRate+1:i*NumCh*SampRate),'uint16');
    end
    fclose(h2);
end

%% digitalin.dat
d = dir('digitalin.dat');
if ~isempty(d)
    disp('Writing digitalin file')
    if RenameOriginalsAsOrig
        inname = 'digitalin.dat_orig';
        outname = 'digitalin.dat';
        movefile(outname,inname)
    else
        inname = 'digitalin.dat';
        outname = 'digitalin_fragment.dat';
    end
    disp('Writing digitalin file')
    NumCh = length(board_dig_in_channels);
    m = memmapfile(inname,'Format','uint16');
    h3 = fopen(outname,'W');
    for i = timeperiod(1)+1:timeperiod(2)
        fwrite(h3,m.Data((i-1)*NumCh*SampRate+1:i*NumCh*SampRate),'uint16');
    end
    fclose(h3);
end

%% time.dat
disp('Writing time file')
if RenameOriginalsAsOrig
    inname = 'time.dat_orig';
    outname = 'time.dat';
    movefile(outname,inname)
else
    inname = 'time.dat';
    outname = 'time_fragment.dat';
end
disp('Writing time file')
m = memmapfile(inname,'Format','int32','writable',false);
h4 = fopen(outname,'W');
fwrite(h4,m.Data(timeperiod(1)*SampRate+1:timeperiod(end)*SampRate),'int32');
fclose(h4);

%% supply.dat
disp('Writing analogin file')
if RenameOriginalsAsOrig
    inname = 'supply.dat_orig';
    outname = 'supply.dat';
    movefile(outname,inname)
else
    inname = 'supply.dat';
    outname = 'supply.dat_fragment.dat';
end
disp('Writing supply file')
NumCh = length(length(supply_voltage_channels));
m = memmapfile(inname,'Format','uint16','writable',false);
h5 = fopen(outname,'W');
fwrite(h5,m.Data(timeperiod(1)*NumCh*SampRate+1:timeperiod(end)*SampRate),'uint16');
fclose(h5);

