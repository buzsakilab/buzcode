function bz_LFPfromDat(basepath,varargin)

% perform lowpass (2 X output Fs) sinc filter on wideband data
% subsample the filtered data and save as a new flat binary
% basename must have basename.dat and basename.xml
% basepath is the full path for basename.dat
% varargin{1} = optional sampling rate, otherwise from XML, otherwise 1KHz
% depends upon iosr tool box https://github.com/IoSR-Surrey/MatlabToolbox
% note that sincFilter was altered to accomodate GPU filtering

%% Contstnts/defaults
default_outFs = 1250;
chunksize = 1e5; % depends on the system... could be bigger I guess

%% Input handling
if ~exist('basepath','var')
basepath = cd;
end
basename = bz_BasenameFromBasepath(basepath);

import iosr.dsp.*

useGPU = false;
if gpuDeviceCount>0
useGPU = true;
end


%get xml
fxml = fullfile(basepath, [basename '.xml']);
if ~exist(fxml,'file')
error('Dat file and/or Xml file does not exist')
end
sizeInBytes = 2; %

syst = bz_getSessionInfo(fxml,'noPrompts',true);
fInfo = dir(fullfile(basepath, [basename '.dat']));
inFs = syst.rates.wideband;
nbChan = syst.nChannels;


%set output sampling rate
if isempty(varargin) && isfield(syst.rates,'lfp') && ~isempty(syst.rates.lfp)
outFs =syst.lfpSampleRate;
elseif ~isempty(varargin)
outFs = varargin{1};

if syst.lfpSampleRate~=outFs
warning(['XML does not match user LFP rate, using user input LFP sampling of: ' num2str(outFs) 'Hz']);
end
else
outFs = default_outFs;
end




ratio =outFs/(inFs/2) /2;
sampleRatio = (inFs/outFs);
ntbuff = find(sampleRatio.*[1:100]>525,1)*sampleRatio;%default filter length  = 1025

%f
nBytes = fInfo.bytes;
nbChunks = floor(nBytes/(nbChan*sizeInBytes*chunksize));



fidI = fopen(fullfile(basepath,[basename,'.dat']), 'r');



fprintf('Extraction of LFP begun \n')

if exist(fullfile(basepath,[basename,'.lfp']))
delete(fullfile(basepath,[basename,'.lfp']))
end

fidout =  fopen(fullfile(basepath,[basename,'.lfp']), 'a');


for ibatch = 1:nbChunks

if mod(ibatch,10)==0
if ibatch~=10
fprintf(repmat('\b',[1 length([num2str(round(100*(ibatch-10)/nbChunks)), ' percent complete'])]))
end
fprintf('%d percent complete', round(100*ibatch/nbChunks));
end

h=waitbar(ibatch/(nbChunks+1));

if ibatch>1
fseek(fidI,((ibatch-1)*(nbChan*sizeInBytes*chunksize))-(nbChan*sizeInBytes*ntbuff),'bof');
dat = fread(fidI,nbChan*(chunksize+2*ntbuff),'int16');
dat = reshape(dat,[nbChan (chunksize+2*ntbuff)]);
else
dat = fread(fidI,nbChan*(chunksize+ntbuff),'int16');
dat = reshape(dat,[nbChan (chunksize+ntbuff)]);
end




DATA = nan(size(dat,1),chunksize/sampleRatio);
for ii = 1:size(dat,1)

d = double(dat(ii,:));

if useGPU
d = gpuArray(d);
end

tmp=  iosr.dsp.sincFilter(d,ratio);
if useGPU

if ibatch==1
DATA(ii,:) = gather_try(int16(real( tmp(sampleRatio:sampleRatio:end-ntbuff))));
else
DATA(ii,:) = gather_try(int16(real( tmp(ntbuff+sampleRatio:sampleRatio:end-ntbuff))));
end

else
if ibatch==1
DATA(ii,:) = int16(real( tmp(sampleRatio:sampleRatio:end-ntbuff)));
else
DATA(ii,:) = int16(real( tmp(ntbuff+sampleRatio:sampleRatio:end-ntbuff)));
end

end

end

fwrite(fidout,DATA(:),'int16');

end



remainder = nBytes/(sizeInBytes*nbChan) - nbChunks*chunksize;
if ~isempty(remainder)

fseek(fidI,((ibatch-1)*(nbChan*sizeInBytes*chunksize))-(nbChan*sizeInBytes*ntbuff),'bof');
dat = fread(fidI,nbChan*(remainder+ntbuff),'int16');
dat = reshape(dat,[nbChan (remainder+ntbuff)]);
DATA = nan(size(dat,1),floor(remainder/sampleRatio));
for ii = 1:size(dat,1)
tmp=  iosr.dsp.sincFilter(double(dat(ii,:)),ratio);

if useGPU

DATA(ii,:) = gather_try(int16(real( tmp(ntbuff+sampleRatio:sampleRatio:end))));
else
DATA(ii,:) = int16(real( tmp(ntbuff+sampleRatio:sampleRatio:end)));
end
end

fwrite(fidout,DATA(:),'int16');


end
close(h);

fclose(fidI);
fclose(fidout);
end

