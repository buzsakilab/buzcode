function keydata = readSvnKeyword(filename,keyword,crop)
%READSVNKEYWORD Read data from a file tagged with an SVN keyword
%   
%   The function takes an input filename and searches that file for a
%   Subversion keyword. The data associated with the keyword are returned
%   in a character array. Keyword data are placed in files automatically by
%   Subversion at each commit using the 'keyword' function. The data can be
%   useful in maintaining an audit trail, or establishing the version used
%   to generate a particular set of results.
%   
%   KEYDATA = IOSR.SVN.READSVNKEYWORD(FILENAME,KEYWORD) returns a string
%   KEYDATA containing data associated with the SVN keyword KEYWORD in a
%   file specified by FILENAME.
%   
%   FILENAME and KEYWORD must be strings specifying a single file and
%   keyword respectively. To read multiple files or use multiple keywords,
%   use BUILD_SVN_PROFILE instead. The function returns an empty string if
%   the keyword is not found.
% 
%   KEYDATA = IOSR.SVN.READSVNKEYWORD(FILENAME,KEYWORD,CROP), with CROP =
%   true (default is false), removes the keyword, and leading and trailing
%   spaces, from the returned string.
%   
%   See also IOSR.SVN.BUILDSVNPROFILE.

%   Copyright 2016 University of Surrey.


    assert(ischar(filename) & ischar(keyword), 'iosr:readSvnKeyword:invalidInputs', 'FILENAME and KEYWORD must be char arrays')

    if nargin < 3
        crop = false;
    end

    keydata = '';

    fid = fopen(filename); % open file
    assert(fid~=-1, 'iosr:readSvnKeyword:invalidFile', ['read_svn_keyword: ''' filename ''' not found'])

    tline = fgetl(fid); % read first line
    while ischar(tline)
        k1 = strfind(tline,['$' keyword ':']); % find keyword
        if ~isempty(k1) % if it is found
            k2 = strfind(tline,'$'); % find the end of the keyword data
            k2 = k2(k2>k1);
            keydata = tline(k1+1:k2-2); % extract the data from the line
            tline = -1; % set tline to numeric
        else
            tline = fgetl(fid); % read next line
        end
    end

    if crop
        keydata = keydata(find(isspace(keydata),1,'first')+1:end);
    end

    fclose(fid); % close file

end
