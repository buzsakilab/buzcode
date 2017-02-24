function fns = FindFiles(globfn, varargin)

% fns = FindFiles(globfn, parameters)
%
% Finds all files that match a wildcard input globfn.
% Based on matlab's dir function.
% Searches all directories under the current directory.
%
% INPUTS
%      globfn -- filename to search for (you can use '*',
%                but not '?', don't use directory names.
%
% OUTPUTS
%       fns -- cell array of files found
%
% PARAMETERS
%       StartingDirectory -- default '.'
%       CheckSubdirs (1/0) -- default 1
%
% ADR 1998
% version L4.2
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

%-----------------
StartingDirectory = '.';
CheckSubdirs = 1;
extract_varargin;
pushdir(StartingDirectory);

%-----------------
thisdirfiles = dir(globfn);
fns = cell(length(thisdirfiles),1);
for iF = 1:length(thisdirfiles)
   fns{iF} = fullfile(pwd,thisdirfiles(iF).name);
end

if CheckSubdirs
   subdirs = dir;
   for iD = 1:length(dir)
      if subdirs(iD).isdir & ...
            strcmp(subdirs(iD).name,'.')==0 & ...       % is not .
            strcmp(subdirs(iD).name,'..')==0            % is not ..      
         pushdir;
         cd(subdirs(iD).name)
         % disp(['FindFiles: searching "',subdirs(iD).name,'".']);
         subfns = FindFiles(globfn);
         popdir;
         fns = [fns; subfns];
      end
   end
end

popdir;


