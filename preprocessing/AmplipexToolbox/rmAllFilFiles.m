function rmAllFilFiles(basedir)

% function rmAllFilFiles(basedir)
%
% This function starts in basedir (or the current directory, if no argument
% is given) and recursively searches through every subdirectory and deletes
% all of the .fil files. The .fil files can be deleted after the ndm_start
% scripts are done.
%
% Luke Sjulson, 2015-05-21


if nargin == 0
    basedir = pwd;
end

% script to remove all .FIL files in the current directory and all
% subdirectories
cd(basedir);
basedir = pwd; % in case the user entered a partial path
k = dir;
kidx = 1;

% ridiculous method to remove names starting with '.'
% want this to run in windows someday if necessary
for idx = 1:length(k)
    if k(idx).name(1)~='.'
        k2(kidx) = k(idx);
        kidx = kidx + 1;
    end
end

if exist('k2')
    if sum([k2.isdir]) > 0 % meaning that it's not in a terminal directory
        whichidxs = find([k2.isdir]);
        for idx = 1:length(whichidxs) % looping over all subdirectories
            rmAllFilFiles([basedir '/' k2(whichidxs(idx)).name]); % recursion
        end
    end
end

cd(basedir);
% fprintf('%s\n', basedir);

if ~isempty(strfind(upper(computer), 'WIN'))
    eval('!del *.fil');
else
    eval('!rm *.fil > /dev/null 2>&1');
end

