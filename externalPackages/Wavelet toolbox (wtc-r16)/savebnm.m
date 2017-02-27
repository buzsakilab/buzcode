function varargout=savebnm(fname,varargin);
% a compact binary format for saving a matrix. (using float)
%
%
% savebnm(fname,M1[,M2,M3,...]);
%
% Can save multidimensional matrices containing: Nans,Imaginary numbers and Inf
%
% see loadbnm
%
% Aslak Grinsted 2004 april
%


% -------------------------------------------------------------------------
%   Copyright (C) 2002-2004, Aslak Grinsted
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.

% TODO: save inputnames also


if (length(strfind(fname,'.'))==0)
    fname=[fname '.bnm'];
end

for ii=1:length(varargin) %validate arguments
    if ~isnumeric(varargin{ii})
        error('BNM format only supports numeric data!');
    end
end

[fid,message]=fopen(fname,'w','ieee-le');
if (fid<0)
    message=strrep(message,'Sorry. No help in figuring out the problem . . .','');
    error(['Could not open the specified file for writing. ' message ' file=' fname])
end

fwrite(fid,'BNM Matrix File');
version=2;

fwrite(fid,[version length(varargin)],'int16');

for ii=1:length(varargin)
    %[ii ftell(fid)]
    %M=varargin{ii};
    sz=size(varargin{ii});
    fwrite(fid,[length(sz) isreal(varargin{ii})],'int16');
    fwrite(fid,sz,'uint32');
    fwrite(fid,varargin{ii},'float'); 
    if ~isreal(varargin{ii})
        fwrite(fid,imag(varargin{ii}),'float'); 
    end
end
fclose(fid);



function savestring(fid,s);
if ~ischar(s)|(size(s,1)~=1)
    error('not a string')
end
fwrite(fid,length(s),'uint32');
fwrite(fid,s);

