function varargout=loadbnm(fname)
% a compact binary matrix format.
%
% [M1,M2,M3,...]=loadbnm(fname);
%
% see savebnm
%
% Aslak Grinsted April 2004


% -------------------------------------------------------------------------
%   Copyright (C) 2002-2004, Aslak Grinsted
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.


if (length(strfind(fname,'.'))==0)
    fname=[fname '.bnm'];
end

skip=(nargout==0); % should it skip loading data and just show info for the file?

header='BNM Matrix File';

[fid,message]=fopen(fname,'r','ieee-le');
if (fid<0)
    message=strrep(message,'Sorry. No help in figuring out the problem . . .','');
    error(['Could not open the specified file. ' message ' file=' fname])
end

hdr=[fread(fid,length(header),'uchar')' ''];
if ~strcmp(hdr,header) 
    error([fname 'is not a ' header ' file!']);
end

version=fread(fid,1,'int16');
if version>2
    warning('File format version > 2')
end
if (version==1)
    nout=1;
else
    nout=fread(fid,1,'int16');
end


nn=min(nargout,nout);
if skip
    disp(sprintf('File: "%s" , Version: %i , MatrixCount : %i\n',fname,version,nout));
    nn=nout;
end

for ii=1:nn
    infos=fread(fid,2,'int16');
    sz=fread(fid,infos(1),'uint32')';
    rl=infos(2); %isreal???
    if skip
        disp(sprintf('Matrix %i,  size: [%s], isimag: %i',ii,num2str(sz),~rl))
        fseek(fid,prod(sz)*(1+(~rl))*4,0);
    else
        varargout{ii}=fread(fid,prod(sz),'float'); 
        if ~rl %if imaginary numbers...
            varargout{ii}=varargout{ii}+fread(fid,prod(sz),'float')*i; 
        end
        varargout{ii}=reshape(varargout{ii},sz);
        %varargout{ii}=M;
    end
end    
    
fclose(fid);






function s=loadstring(fid)
ll=fread(fid,1,'uint32');
s=[fread(fid,ll,'uchar')' ''];
