function bz_LFPfromDat(basepath,varargin)
% perform lowpass (2 X output Fs) sinc filter on wideband data
% subsample the filtered data and save as a new flat binary
% basename must have basename.dat and basename.xml
% basepath is the full path for basename.dat
%
% note that sincFilter was altered to accomodate GPU filtering
%
%INPUTS
%   basePath    path where the recording files are located
%               where basePath is a folder of the form: 
%                   whateverPath/baseName/
%
%               Assumes presence of the following files:
%                   basePath/baseName.dat
%                   -or-
%                   basePath/amplifier.dat
%
%                   (optional parameters files)
%                   basePath/baseName.xml
%                   basePath/baseName.sessionInfo.mat
%
%               If basePath not specified, tries the current working directory.
%
%   (options)
%       'outFs'         (default: 1250) downsampled frequency of the .lfp
%                       output file. if no user input and not specified in
%                       the xml, use default
%       'lopass'        (default: 450) low pass filter frequency 
%       'noPrompts'     (default: true) prevents prompts about
%                       saving/adding metadata
%
%
%OUTPUT
%   Creates file:   basePath/baseName.lfp
%
%   If no sessionInfo.mat file previously exists, creates one with 
%   the information from the .xml file, with the .lfp sampling frequency 
%   and the lowpass filter used.
%
%
%Dependency: iosr tool box https://github.com/IoSR-Surrey/MatlabToolbox
%
%SMckenzie, BWatson, DLevenstein 2018
%% Input handling
if ~exist('basepath','var')
    basepath = pwd;
end
basename = bz_BasenameFromBasepath(basepath);

defaultoutFS = 1250; %used later

p = inputParser;
addParameter(p,'noPrompts',true,@islogical);
addParameter(p,'outFs',[],@isnumeric);
addParameter(p,'lopass',450,@isnumeric);
parse(p,varargin{:})
noPrompts = p.Results.noPrompts;
outFs = p.Results.outFs;
lopass = p.Results.lopass;

import iosr.dsp.*

useGPU = false;
try
    if gpuDeviceCount>0
        useGPU = true;
    end
end
sizeInBytes = 2; %

%% files check
fxml = fullfile(basepath, [basename '.xml']);
fsessioninfo = fullfile(basepath,[basename,'.sessionInfo.mat']);
fdat = fullfile(basepath,[basename,'.dat']);
flfp = fullfile(basepath,[basename,'.lfp']);

%If there's already a .lfp file, make sure the user wants to overwrite it
if exist(flfp,'file')
    overwrite = input([basename,'.lfp already exists. Overwrite? [Y/N] '],'s');
    switch overwrite
        case {'y','Y'}
            delete(flfp)
        case {'n','N'}
            return
        otherwise
            error('Y or N please...')
    end
end

%Check the dat
if ~exist(fdat,'file')
    fdat = fullfile(basepath,'amplifier.dat'); %Try amplifier.dat
    if ~exist(fdat,'file')
        error('Dat file does not exist')
    end
    
end
fInfo = dir(fullfile(basepath, [basename '.dat']));

%Get the metadata
if ~exist(fxml,'file') && ~exist(fsessioninfo,'file')
    warning('No xml or sessionInfo file, using defaults and creating a minimal sessionInfo file')
else
    %Get everything from the xml/sessionInfo
    sessionInfo = bz_getSessionInfo(basepath,'noPrompts',noPrompts);
    inFs = sessionInfo.rates.wideband;
    nbChan = sessionInfo.nChannels;
    
    %set output sampling rate from xml, user input    
    if ~isempty(outFs)          %If user input - priority (keep from above)
        outFs = outFs;          %redundant, clearly.
    elseif isfield(sessionInfo,'lfpSampleRate') %If not, use from xml
        outFs = sessionInfo.lfpSampleRate;    
    else                                        %If not in xml, use default
        outFs = defaultoutFS;
    end
end
sessionInfo.lfpSampleRate = outFs;
sessionInfo.LFPLoPassFreq = lopass;
save(fsessioninfo,'sessionInfo');  %Save the sessioninfo with the parameters used


if lopass> outFs/2
    warning('low pass cutoff beyond Nyquist')
end
 
ratio =lopass/(inFs/2) ;
sampleRatio = (inFs/outFs);

%% Set Chunk and buffer size at even multiple of sampleRatio
chunksize = 1e5; % depends on the system... could be bigger I guess
if mod(chunksize,sampleRatio)~=0
    chunksize = chunksize + sampleRatio-mod(chunksize,sampleRatio);
end

%ntbuff should be even multiple of sampleRatio
ntbuff = 525;  %default filter size in iosr toolbox
if mod(ntbuff,sampleRatio)~=0
    ntbuff = ntbuff + sampleRatio-mod(ntbuff,sampleRatio);
end

nBytes = fInfo.bytes;
nbChunks = floor(nBytes/(nbChan*sizeInBytes*chunksize));

%% GET LFP FROM DAT!
fidI = fopen(fdat, 'r');
fprintf('Extraction of LFP begun \n')
fidout = fopen(flfp, 'a');

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

disp(' ........baseName.lfp file created! Huzzah!')
end

