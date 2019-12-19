function [t,amps,data,aux] = read_intan_data(filename)

% [t,amps,data,aux] = read_intan_data
%
% Opens file selection GUI to select and then read data from an Intan 
% amplifier data file (*.int).
%
% t = time vector (in seconds)
% amps = vector listing active amplifier channels
% data = matrix of electrode-referred amplifier signals (in microvolts)
% aux = matrix of six auxiliary TTL input signals
%
% Example usage:
%  >> [t,amps,data,aux] = read_intan_data;
%  >> plot(t,data(:,1));
%
% Version 1.1, June 26, 2010
% (c) 2010, Intan Technologies, LLC
% For more information, see http://www.intantech.com
% For updates and latest version, see http://www.intantech.com/software.html
%
% 06-22-10 Added GUI file selection and optimized: Craig Patten, Plexon, Inc.

% use MATLAB predefined gui uigetfile to select the file(s) to analyze
if ~exist('filename','var')
    [file, path, filterindex] = uigetfile('*.int','Select a .int file','MultiSelect', 'off');
    filename = [path,file];
end

fid = fopen(filename, 'r');

% Read first three header bytes encoding file version
for i=1:3
    header(i) = fread(fid, 1, 'uint8');
end

if (header(1) ~= 128)
    error('Improper data file format.');
end

if (header(2) ~= 1 || header(3) ~= 1)
    warning('Data file version may not be compatible with this m-file.');
end

% Now see which amplifier channels are saved in this file.
for i=1:64
    amp_on(i) = fread(fid, 1, 'uint8');
end

num_amps = sum(amp_on);

% Create a list of amplifier channels in this file.
amps = zeros(1,num_amps);
index = 1;
for i=1:64
    if (amp_on(i) == 1)
        amps(index) = i;
        index = index + 1;
    end
end

% Now search for the end of the file to find out the length of the data.
% t_count = 0;
% while (~feof(fid))
%    fread(fid, 1+4*num_amps, 'uint8'); 
%    t_count = t_count + 1;
% end
% t_count = t_count - 1;
% t_max = t_count/25000;

%-----------------------------------
% replace above code with a more efficient method CDP 06-24-10
s = dir(filename);
filesize = s.bytes;
t_count = (filesize - 67)/(num_amps*4 + 1);
t_max = t_count/25000;
%-----------------------------------

% print channel (singular) when there is only one channel! CDP 06-24-10
if num_amps == 1;
    fprintf(1, '\nData file contains %0.2f seconds of data from %d amplifier channel.\n', t_max, num_amps);
    fprintf(1, 'Channel: ');
else
    fprintf(1, '\nData file contains %0.2f seconds of data from %d amplifier channels.\n', t_max, num_amps);
    fprintf(1, 'Channels: ');
end

for i=1:num_amps
    fprintf(1, '%d ', amps(i));
end
fprintf(1, '\n\n');

% Pre-allocate large data matrices.
aux = zeros(t_count,6,'uint8');
t = (0:1:(t_count-1))/25000;
t = t';

%--------------------------------------
% Replace code code below with much faster code CDP 06-24-10
% Go back to the beginning of the file...
frewind(fid);

% ...skip the header this time...
fread(fid, 3+64, 'uint8');

% allocate space to read the entire file
data2 = zeros((filesize-67),1,'uint8');
% read the entire file
data2 = fread(fid,(filesize-67),'uint8=>uint8');

% extract the digital data
aux_data = data2((num_amps*4)+1:num_amps*4+1:filesize-67);

% extract individual bits
aux = [bitget(aux_data,6),bitget(aux_data,5),bitget(aux_data,4),bitget(aux_data,3),bitget(aux_data,2),bitget(aux_data,1)];
clear aux_data;

% delete the digital data
data2((num_amps*4)+1:num_amps*4+1:filesize-67) = [];

% convert the remaining data from bytes to single
data2 = typecast(data2,'single');

data = zeros(t_count,num_amps);
% de-mux the channels
for ind = 1:num_amps
    data(:,ind) = data2(ind:num_amps:length(data2));
end
%---------------------------------------

% % Go back to the beginning of the file...
% frewind(fid);
% 
% % ...skip the header this time...
% fread(fid, 3+64, 'uint8');
% 
% % ...and read all the data.
% fprintf(1, 'Reading data...  (This may take a while.)\n\n');
% for i=1:t_count
%     for j=1:num_amps
%         data(i,j) = double(fread(fid, 1, 'float32'));
%     end
%     
%     aux_byte = fread(fid, 1, 'uint8');
%     
%     % Decode auxiliary TTL inputs
%     if aux_byte >= 32
%         aux(i,6) = 1;
%         aux_byte = aux_byte - 32;
%     end
%     if aux_byte >= 16
%         aux(i,5) = 1;
%         aux_byte = aux_byte - 16;
%     end
%     if aux_byte >= 8
%         aux(i,4) = 1;
%         aux_byte = aux_byte - 8;
%     end
%     if aux_byte >= 4
%         aux(i,3) = 1;
%         aux_byte = aux_byte - 4;
%     end
%     if aux_byte >= 2
%         aux(i,2) = 1;
%         aux_byte = aux_byte - 2;
%     end
%     if aux_byte >= 1
%         aux(i,1) = 1;
%         aux_byte = aux_byte - 1;
%     end        
% end

% Close file, and we're done.
fclose(fid);
