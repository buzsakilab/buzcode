function [returnVar,msg] = ReRefFilFile(fname,nbChan,refChan,rerefChan)

% USAGE:
%     ReRefFilFile(fbasename, nbChan, refChan, rerefChan)
%     This function substracts one channel or the median of a series of
%     channels (the 'reference' channel) from other channels. 
%
%     Note: this uses 1-based indexing, not zero-based, so you have to take
%     the channel number in neuroscope and add one to it, e.g. channel 72 
%     in neuroscope is channel 73 here.
%
% INPUTS:
%     fname: fil file name. Could be also a file basename not finishing
%     with '.fil' and the program will be then executed for all the fil
%     files begining with this file basename.
%     nbChan: total number of channels in fil file
%     refChan: reference channel (could be a vector, in this case the median will be the new reference)
%     rerefChan: vector of channels to rereference
% 
% Adrien Peyrache 2013

if ~strcmp(fname(end-2:end),'fil')
    datFiles = dir([fname '*.fil']);
    for ii=1:length(datFiles)
        ReRefFilFile(datFiles(ii).name,nbChan,refChan,rerefChan)
    end

else
    fprintf('ReReferencing %s\n',fname)
    try
        infoFile = dir(fname);

        chunk = 1e6;
        nbChunks = floor(infoFile.bytes/(nbChan*chunk*2));
        warning off
        if nbChunks==0
            chunk = infoFile.bytes/(nbChan*2);
        end

        for ix=0:nbChunks-1
            m = memmapfile(fname,'Format','int16','Offset',ix*chunk*nbChan*2,'Repeat',chunk*nbChan,'writable',true);
            d = m.Data;
            d = reshape(d,[nbChan chunk]);
            ref = d(refChan,:);
            if length(refChan)>1
                ref = median(double(ref));
            end
            ref = repmat(ref,[length(rerefChan) 1]);
            d(rerefChan,:) = d(rerefChan,:)-int16(ref);

            m.Data = d(:);
            clear d m
        end
        %close(h)


        newchunk = infoFile.bytes/(2*nbChan)-nbChunks*chunk;

        if newchunk
            m = memmapfile(fname,'Format','int16','Offset',nbChunks*chunk*nbChan*2,'Repeat',newchunk*nbChan,'writable',true);
            d = m.Data;
            d = reshape(d,[nbChan newchunk]);
            ref = d(refChan,:);
            if length(refChan)>1
                ref = median(ref);
            end
            ref = repmat(ref,[length(rerefChan) 1]);
            d(rerefChan,:) = d(rerefChan,:)-int16(ref);
            m.Data = d(:);
         clear d m
        end
        warning on
        returnVar = 1;
        msg = '';

    catch
        fprintf(['Error occurred in processing ' fname '. File not rereferenced.\n']);
        keyboard
        returnVar = 0;
        msg = lasterr; 
    end
    clear m
end