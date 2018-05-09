function bz_LFPfromDat(basepath,varargin)
% perform lowpass (2 X output Fs) sinc filter on wideband data
% subsample the filtered data and save as a new flat binary
% basename must have basename.dat and basename.xml
% basepath is the full path for basename.dat
% depends upon iosr tool box https://github.com/IoSR-Surrey/MatlabToolbox
% note that sincFilter was altered to accomodate GPU filtering
% 
% option inputs
% 
%   outFs = downsampled frequency (1250)
%   lopass = low pass cut off for sinc filter (400)
%%

%% Input handling
if ~exist('basepath','var')
    basepath = pwd;
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

%Contstnts/defaults
outFs = 1250;
lopass = 450;
chunksize = 1e5; % depends on the system... could be bigger I guess
%chunk should be even multiple of sampleRatio

%set output sampling rate from xml
if isfield(syst.rates,'lfp') && ~isempty(syst.rates.lfp)
    outFs =syst.lfpSampleRate;
end

%Parsing the input options
for i = 1:2:length(varargin)
    switch varargin{i}
        
        case 'outFs'    %User set output frequency
                outFs = varargin{i+1};
                if isfield(syst,'lfpSampleRate') &&  syst.lfpSampleRate~=outFs
                    warning(['XML does not match user LFP rate, using user input LFP sampling of: ' num2str(outFs) 'Hz']);
                end
                
        case 'lopass'   %User set lowpass frequency
                lopass = varargin{i+1};
                    
    end
end

if lopass> outFs/2
    warning('low pass cutoff beyond Nyquist')
end

    
ratio =lopass/(inFs/2) ;
sampleRatio = (inFs/outFs);


if mod(chunksize,sampleRatio)~=0
chunksize = chunksize + sampleRatio-mod(chunksize,sampleRatio);
end


%ntbuff should be even multiple of sampleRatio
ntbuff = 525;  %default filter size in iosr toolbox

if mod(ntbuff,sampleRatio)~=0
ntbuff = ntbuff + sampleRatio-mod(ntbuff,sampleRatio);
end



%f
nBytes = fInfo.bytes;
nbChunks = floor(nBytes/(nbChan*sizeInBytes*chunksize));


%If there's already a .lfp file, make sure the user wants to overwrite it
if exist(fullfile(basepath,[basename,'.lfp']),'file')
    overwrite = input([basename,'.lfp already exists. Overwrite? [Y/N]']);
    switch overwrite
        case {'y','Y'}
            delete(fullfile(basepath,[basename,'.lfp']))
        case {'n','N'}
            return
        otherwise
            error('Y or N please...')
    end
end

fidI = fopen(fullfile(basepath,[basename,'.dat']), 'r');
fprintf('Extraction of LFP begun \n')
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
        d = double(dat(ii,:));
        if useGPU
            d = gpuArray(d);
        end
        
        tmp=  iosr.dsp.sincFilter(d,ratio);
        
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

