%function [data OrigIndex]= LoadBinary_bw(FileName, Channels, nChannels, varargin)
%   Channels - list of channels to load starting from 1
%   nChannels - number of channels in a file, will be read from par/xml file
% if present
%   intype/outtype - data types in the file and in the matrix to load to
% by default assume input file is eeg/dat = int16 type (short int), and
% output is single (to save space) unless it is version 6.5 that cannot
% handle single type
%   method: (1,2,  3 or 4) differes by the way the loading is done - just for
% efficiency purposes some are better then others, default =2;
% method 2 works with buffers-works even for huge files. other methods
% don't work so far ..
% NB: method =3 allows to load data from within certain time epochs , give
% in variable Periods : [beg1 end1; beg2 end2....] (in sampling rate of the
% file you are loading, so if you are loading eeg - then Periods should be
% in eeg samples
% OrigIndex then returns the original samples index that samples in Data correspond
% to , so that you can use it for future spikes and other point process
% analysis
% NB: for method 4 for efficiency and historical reasons output is nCh x nT 
% complaints to : Anton


function [data OrigIndex]= LoadBinary_bw(FileName, Channels, varargin)
if ~FileExists(FileName)
    error('File %s does not exist or cannot be open\n',FileName);
end

lastdot =strfind(FileName,'.');
FileBase=FileName(1:lastdot(end)-1);
if FileExists([FileBase '.xml']) || FileExists([FileBase '.par'])
    Par = LoadPar([FileBase]);
    nChannels = Par.nChannels;
else
    nChannels = 1;
end

[nChannels, method, intype, outtype,Periods,Resample] = DefaultArgs(varargin,{nChannels,2,'int16','double',[],1});

if ~nChannels error('nChannels is not specified and par/xml file is not present'); end

ver = version; ver = str2num(ver(1));
if ver<7 outtype ='double';end

PrecString =[intype '=>' outtype];
fileinfo = dir(FileName);
% get file size, and calculate the number of samples per channel
nSamples = ceil(fileinfo(1).bytes /datatypesize(intype) / nChannels);

if method<2 & ~isempty(Periods)
    error('this method does not perform (yet) selective loading with Periods, use method 3 or 4. Bug me to implement it ! :)');
end
if method==3 & isempty(Periods)
    method=2;
    fprintf('method 3 is replaced by method 2 which uses buffering');
end

%have not fixed the case of method 1 for periods extraction
if method<5

    filept = fopen(FileName,'r');

    if ~isempty(Periods)
        %        method=3;
        if Resample>1
            data = feval( outtype, zeros( length(Channels), sum( ceil((diff(Periods,1,2)+1)/nChannels/Resample) ) ) );
        else
            data = feval( outtype, zeros( length(Channels), sum(diff(Periods,1,2)+1) ) );
        end
    else
        Periods = [1 nSamples];
        if Resample>1
            data = feval(outtype,zeros(length(Channels), ceil(nSamples/Resample)));
        else
            data = feval(outtype,zeros(length(Channels), nSamples));
        end
    end
end

switch method
    case 1

        %compute continuous patches in chselect
        %lazy to do circular continuity search - maybe have patch [nch 1 2 3]
        %sort in case somebody didn't
        [Channels ChanOrd]= sort(Channels(:)');
        dch = diff([Channels(1) Channels]);
        PatchInd = cumsum(dch>1)+1;
        PatchLen = hist(PatchInd,unique(PatchInd));
        PatchSkip = (nChannels - PatchLen)*datatypesize(intype);
        nPatch = length(unique(PatchInd));

        for ii=1:nPatch
            patchch = find(PatchInd==ii);
            patchbeg = Channels(patchch(1));6
            PatchPrecString = [num2str(PatchLen(ii)) '*' PrecString];
            fseek(filept,(patchbeg-1)*datatypesize(intype),'bof');
            data(patchch,:) = fread(filept,[PatchLen(ii) nSamples],PatchPrecString,PatchSkip(ii));

        end
        % put them back in the order they were in Channels argument
        data = data(ChanOrd,:);

        
    case 2 %old way - buffered, now deals with periods as well
        OrigIndex = [];
        nPeriods = size(Periods,1);
        buffersize = 400000;
        if Resample>1 buffersize = round(buffersize/Resample)*Resample; end
        totel=0;
        for ii=1:nPeriods
            numel=0;
            numelm=0;
            Position = (Periods(ii,1)-1)*nChannels*datatypesize(intype);
            ReadSamples = diff(Periods(ii,:))+1;
            fseek(filept, Position, 'bof');
            while numel<ReadSamples 
                if numel==ReadSamples break; end
                [tmp,count] = fread(filept,[nChannels,min(buffersize,ReadSamples-numel)],PrecString);
                data(:,totel+1:totel+ceil(count/nChannels/Resample)) = tmp(Channels,1:Resample:end);
                clear tmp
               
                numel = numel+count/nChannels;
                totel = totel+ceil(count/nChannels/Resample);
            end
            
            OrigIndex = [OrigIndex; Periods(ii,1)+[0:Resample:ReadSamples-1]'];
        end

    case 3 % for full periods extraction, not buffered, use method 2 if periods are large.
        % OBSOLETE!!!
        nPeriods = size(Periods,1);
        numel=0;
        OrigIndex = [];
        for ii=1:nPeriods
            Position = (Periods(ii,1)-1)*nChannels*datatypesize(intype);
            ReadSamples = diff(Periods(ii,:))+1;
            fseek(filept, Position, 'bof');
            
            [tmp count]= fread(filept, [nChannels, ReadSamples], PrecString);
            if count/nChannels ~= ReadSamples
                error('something went wrong!');
            end
            if Resample>1
                numelm = ceil(count/nChannels/Resample);
            else
                numelm = count/nChannels;
            end
            data(:,numel+1:numel+numelm) = tmp(Channels,1:Resample:end);
            numel = numel+count/nChannels;
            OrigIndex = [OrigIndex; Periods(ii,1)+[0:Resample:ReadSamples-1]'];
        end
        
    case 4 %new way - with memmapfile

        if isempty(Periods)
            mmap = memmapfile(FileName, 'format',{intype [nChannels nSamples] 'x'},'offset',0,'repeat',1);
            data = mmap.Data.x(Channels,1:Resample:end);
        else
            %data = feval(outtype,zeros(length(Channels), sum(diff(Periods,1,2)))+size(Periods,1));
            nPeriods = size(Periods,1);
            data = [];
            OrigIndex = [];
            for ii=1:nPeriods
                Position = (Periods(ii,1)-1)*nChannels*datatypesize(intype);
                ReadSamples = diff(Periods(ii,:))+1;
                mmap = memmapfile(FileName, 'format',{intype [nChannels ReadSamples] 'x'},'offset',Position,'repeat',1);
%                data = [data mmap.Data.x(Channels,Periods(ii,1):Periods(ii,2))];
                data = [data mmap.Data.x(Channels,1:Resample:end)];
                OrigIndex = [OrigIndex; Periods(ii,1)+[0:Resample:ReadSamples-1]'];
            end
        end
        data = cast(data,outtype);
end
if method<4
    fclose(filept);
end


