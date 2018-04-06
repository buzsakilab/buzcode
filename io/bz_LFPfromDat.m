
function bz_LFPfromDat(basename,basepath,varargin)

% perform lowpass (2 X output Fs) sinc filter on wideband data
% subsample the filtered data and save as a new flat binary
% basename must have basename.dat and basename.xml
% basepath is the full path for basename.dat
% varargin{1} = optional sampling rate, otherwise from XML, otherwise 1KHz
% depends upon iosr tool box https://github.com/IoSR-Surrey/MatlabToolbox
% note that sincFilter was altered to accomodate GPU filtering



import iosr.dsp.*

useGPU = false;
if gpuDeviceCount>0
    useGPU = true;
end


%get xml
fxml = [basepath '/' basename '.xml'];
if ~exist(fxml,'file')
    error('Dat file and/or Xml file does not exist')
end
sizeInBytes = 2; %
chunk = 1e5; % depends on the system... could be bigger I guess

syst = LoadXml(fxml);
fInfo = dir([basepath '/' basename '.dat']);
inFs = syst.SampleRate;
nbChan = syst.nChannels;


%set output sampling rate
if isempty(varargin) && isfield(syst.lfpSampleRate) && ~isempty(syst.lfpSampleRate)
    outFs =syst.lfpSampleRate;
elseif ~isempty(varargin)
    outFs = varargin{1};
    
    if syst.lfpSampleRate~=outFs
        warning(['XML does not match user LFP rate, using user input LFP sampling of: ' num2str(outFs) 'Hz']);
    end
else
    outFs = 1000;
end




ratio =outFs/(inFs/2) /2;
sampleRatio = (inFs/outFs);
ntbuff = find(sampleRatio.*[1:100]>525,1)*sampleRatio;%default filter length  = 1025

%f
nBytes = fInfo.bytes;
nbChunks = floor(nBytes/(nbChan*sizeInBytes*chunk));



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
        fseek(fidI,((ibatch-1)*(nbChan*sizeInBytes*chunk))-(nbChan*sizeInBytes*ntbuff),'bof');
        dat = fread(fidI,nbChan*(chunk+2*ntbuff),'int16');
        dat = reshape(dat,[nbChan (chunk+2*ntbuff)]);
    else
        dat = fread(fidI,nbChan*(chunk+ntbuff),'int16');
        dat = reshape(dat,[nbChan (chunk+ntbuff)]);
    end
    
    
    
    
    DATA = nan(size(dat,1),chunk/sampleRatio);
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



remainder = nBytes/(sizeInBytes*nbChan) - nbChunks*chunk;
if ~isempty(remainder)
    
    fseek(fidI,((ibatch-1)*(nbChan*sizeInBytes*chunk))-(nbChan*sizeInBytes*ntbuff),'bof');
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
