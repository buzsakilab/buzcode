function [rev,file] = headRev(folders,strs,rec)
%HEADREV Retrieve the head revision for specified files
% 
%   This function finds the most recently committed file in the specified
%   folders(s) and returns its revision and filename. The files must
%   contain the 'Revision' and 'Date' keywords.
% 
%   REV = IOSR.SVN.HEADREV(FOLDERS) returns the head revision (i.e. the
%   highest revision) for all files in FOLDERS.
%
%   REV = IOSR.SVN.HEADREV(FOLDERS,STRS) allows additional filter strings
%   STRS to be specified. See GETCONTENTS for details of permitted filter
%   strings.
% 
%   REV = IOSR.SVN.HEADREV(FOLDERS,STRS,REC), with REC = true, allows
%   folders to be searched recursively (default is false).
%
%   [REV,FILE] = IOSR.SVN.HEADREV(...) returns the filename to FILE of the
%   file with the highest revision.
%
%   FOLDERS and STRS may be a string specifying a single occurrence, or a
%   cell array of strings specifying multiple occurrences. This function
%   requires that the 'Revision' keyword is used in the searched functions.
% 
%   Note that if the folders include files from an externally-defined
%   repository (which has been updated more recently than the native
%   repository), a misleading revision number may be presented.
% 
%   See also IOSR.GENERAL.GETCONTENTS, IOSR.SVN.BUILDSVNPROFILE,
%            IOSR.SVN.READSVNKEYWORD.

%   Copyright 2016 University of Surrey.

    keyword = 'Date'; % use this data to sort files (doubtful that any others would work, without considerable parsing)

    if ischar(folders)
        folders = cellstr(folders);
    end

    if nargin > 1
        if ischar(strs)
            strs = cellstr(strs);
        end
        if isempty(strs)
            strs = {'files'};
        end
    else
        strs = {'files'};
    end

    % add subfolders if recursion is requested
    if nargin > 2
        if rec
            for f = folders
                folders = [folders; iosr.general.getContents(char(f),'filter','folders','rec',true,'path','full')]; %#ok<AGROW>
            end
        end
    end

    % get revision info for files
    svn_profile = iosr.svn.buildSvnProfile(folders,keyword,strs);

    % remove the keyword and convert revision strings to numbers
    svn_profile(:,2) = cellfun(@(x) (x(length(keyword)+3:end)),svn_profile(:,2),'uni',false);

    % remove NaN (files that don't have the revision keyword)
    IX1 = cellfun(@isempty,svn_profile(:,2));
    svn_profile = svn_profile(~IX1,:);

    % sort by revision number
    [~,IX2] = sort(svn_profile(:,2));
    svn_profile = svn_profile(IX2,:);

    % get highest revision number
    rev = iosr.svn.readSvnKeyword(svn_profile{end,1},'Revision');

    % return corresponding filename, if requested
    if nargout>1
        file = svn_profile{end,1};
    end

end
