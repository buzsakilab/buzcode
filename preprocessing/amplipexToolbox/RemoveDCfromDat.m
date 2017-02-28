function [returnVar,msg] = RemoveDCfromDat(fname,nbChan)

% USAGE:
%     RemoveDCfromDat(fbasename,nbChan)
%     This function removes DC from dat files by computing the average of
%     the first 1e6 samples (or less if file is smaller)
% INPUTS:
%     fname: dat file name
%     nbChan: total number of channels in dat file
% 
% Adrien Peyrache 2011

fprintf('Removing baseline from %s\n',fname)
try
    infoFile = dir(fname);
    
    chunk = 1e6;
    nbChunks = floor(infoFile.bytes/(nbChan*chunk*2));
    warning off
    if nbChunks==0
        chunk = infoFile.bytes/(nbChan*2);
    end
    m = memmapfile(fname,'Format','int16','Repeat',chunk*nbChan,'writable',true);
    d = m.Data;
    d = reshape(d,[nbChan chunk]);
    meanD = mean(d,2);
    d = d-int16(meanD*ones(1,chunk));
    m.Data = d(:);
    clear d m
    
    for ix=1:nbChunks-1
    %    h=waitbar(ix/(nbChunks-1));
        m = memmapfile(fname,'Format','int16','Offset',ix*chunk*nbChan*2,'Repeat',chunk*nbChan,'writable',true);
        d = m.Data;
        d = reshape(d,[nbChan chunk]);
        d = d-int16(meanD*ones(1,chunk));
        m.Data = d(:);
        clear d m
    end
    %close(h)
    
    
    newchunk = infoFile.bytes/(2*nbChan)-nbChunks*chunk;
   
    if newchunk
        m = memmapfile(fname,'Format','int16','Offset',nbChunks*chunk*nbChan*2,'Repeat',newchunk*nbChan,'writable',true);
        d = m.Data;
        d = reshape(d,[nbChan newchunk]);
        d = d-int16(meanD*ones(1,newchunk));
        m.Data = d(:);
     clear d m
    end
    warning on
    returnVar = 1;
    msg = '';
    
catch
%     keyboard
    returnVar = 0;
    msg = lasterr; 
end
clear m
