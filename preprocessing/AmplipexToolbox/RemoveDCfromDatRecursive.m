function RemoveDCfromDatRecursive(basedir, fakeRunYN)

% function RemoveDCfromDatRecursive(basedir, fakeRunYN)
%
% This function recursively searches basedir and all subdirectories and
% finds all of the .dat files. The number of channels is read in
% automatically from the corresponding .meta files, and RemoveDCfromDat is
% run on each .dat file. When the DC shift is removed from whatever.dat, 
% a corresponding whatever.DC file is added to that directory to mark the 
% .dat file as completed. When RemoveDCfromDatRecursive is run a second 
% time on the same directory tree, any .dat file with a corresponding .DC 
% file is skipped. 
%
% If you have .dat files that have already had the DC shift subtracted, and
% you want to mark them as done by making .DC files, set fakeRunYN to 'y'.
%
% One other caveat: this only subtracts DC shift from .dat files that have
% corresponding .meta files, so if you have .dat files that were not
% directly recorded by amplirec (e.g. concatenated .dat files), the DC
% shift will not be subtracted.
%
% Example:
%
% >> cd whereYourDATfilesAre
% >> RemoveDCfromDat
%
%  -or-
% 
% >> RemoveDCfromDat('whereYourDATfilesAre')
% 
%  -or-
%
% >> RemoveDCfromDat('whereYourDATfilesAreWithDCalreadySubtracted', 'y') 
%     (marks .dat files as completed without subtracting DC)
%
% Luke Sjulson, 2013-12-11



% clear all
% close all
% clc
% 
% 
% % for testing
% 
% %cd /data/M779/2013-12-02/site1;
% %basename = 'M779_n01';
% basedir = '/data/M779/2013-11-27';

% actual function
if nargin<1
    basedir = '.';
end

if nargin<2
    fakeRunYN = 'n';
end

cd(basedir);
basedir = pwd; % this is to convert a relative path into an absolute path
basedir(basedir=='\') = '/'; % because of windows's stupid backslashes
dirlist = dir;
isdir = find([dirlist.isdir]);

% loop over all subdirectories of basedir recursively
if length(isdir)>2
    for diridx = 3:length(isdir)
        newbasedir = sprintf('%s/%s', basedir, dirlist(isdir(diridx)).name);
        RemoveDCfromDatRecursive(newbasedir, fakeRunYN);
        cd(basedir);
    end
end

% now look for dat files
datlist = dir('*.dat');
if ~isempty(datlist)
    for datidx = 1:length(datlist)
        fname = datlist(datidx).name;
        basename = fname(1:end-4);
        
        % check if the DC has already been removed
        DClist = dir(sprintf('%s.DC', basename));
        if isempty(DClist)
            
            % read in number of channels from .meta file
            fid = fopen(sprintf('%s.meta', basename), 'rt');
            if fid==-1
                warning(sprintf('Could not open .meta file for %s/%s.dat ...skipping', basedir, basename));
                
            else
                textstuff = textscan(fid, '%s%s', 'delimiter', '=');
                nChannelsTxt = textstuff{2}(strcmp('Number of recorded channels ', textstuff{1}));
                nChannels = str2num(nChannelsTxt{1});
                fclose(fid);
                
                % run RemoveDCfromDat
                fprintf('%s:\n', basedir);
                
                
                if upper(fakeRunYN) ~= 'Y'
                    fprintf('RemoveDCfromDat(''%s'', %d);\n\n', fname, nChannels);
                    eval(sprintf('RemoveDCfromDat(''%s'', %d);\n\n', fname, nChannels));
                else
                    fprintf('Fake RemoveDCfromDat(''%s'', %d) - .DC file written;\n\n', fname, nChannels);
                end
                
                % write DC file
                DCfid = fopen(sprintf('%s.DC', basename), 'at');
                % DCfid = 1; % write to screen
                fprintf(DCfid, 'RemoveDCfromDatRecursive just puts this file here so you know that the DC has been subtracted.\n');
                fclose(DCfid);
            end
            
            
        else
            fprintf('%s:\n', basedir);
            fprintf('%s DC file detected\n\n', fname);
        end
    end
else
    disp([pwd ':']);
    fprintf('no dat files here.\n\n');
end


