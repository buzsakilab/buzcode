function logAnalysisFile(AnalysisFileName, writePath)
% The purpose of this utility is to provide a record of analyses run by
% writing out to a log file the state of the code that was used to analyze
% the data. It writes the log file into the current directory if only a
% filename is provided. An absolute path is required to write out to an
% arbitrary location.
%
% This utility is designed to be placed within an .m file used for
% analysis.
%
% dependencies: mfilename
%
% author: John D. Long II, PhD   contact: jlong29@gmail.com
%
%%%%%%%%%%%%%%
%%% INPUTS %%%
%%%%%%%%%%%%%%
% AnalysisFileName: the output of mfilename('fullname') called from within
%   the .m file
% writePath (optional): a user specified path for writing the log file.
%
%%%%%%%%%%%%%%%
%%% OUTPUTS %%%
%%%%%%%%%%%%%%%
% A .log file entitled mfilename-date.log written to the location 
%
% Example usage:
%   logAnalysisFile(mfilename('fullpath') or 
%   logAnalysisFile(mfilename('fullpath','~/mnt/Analysis/')

% input check
if nargin < 2 || isempty(writePath)
    writePath = [];
end

% assumes output from mfilename
info         = dir([AnalysisFileName '.m']);
[~,filename] = fileparts([AnalysisFileName '.m']);

% Check for log file and create or update as required
if ~isempty(writePath)
    % check if writePath is a valid directory
    if exist(writePath,'dir')
        % check if writePath ends in filesep
        if ~strcmp(writePath(end),filesep)
            writePath = [writePath filesep];
        end
        % check for pre-existing log file in writePath
        loginfo = dir(sprintf('%s%s.log',writePath,filename));
        if isempty(loginfo)
            fprintf(1,'First Use Of Analysis File: Writing Log.\n');
            fid1 = fopen(sprintf('%s%s.log',writePath,filename),'w+');
        else
            % let's check if we need to update this file
            fid1  = fopen(sprintf('%s%s.log',writePath,filename),'r');
            tline = fgetl(fid1);
            if strcmp(tline,info.date)
                fprintf(1,'Analysis File Up To Date.\n');
                return
            else
                fprintf(1,'Analysis File Changed: Writing Log.\n');
                fclose(fid1);
                fid1 = fopen(sprintf('%s%s.log',writePath,filename),'w+');
            end
        end
    else
        
        fprintf(1,'InputWarning: Log File Written out to launch directory.\n');
        
        % check for pre-existing log file current directory
        loginfo = dir(sprintf('%s.log',filename));
        if isempty(loginfo)
            fprintf(1,'First Use Of Analysis File: Writing Log.\n');
            fid1 = fopen(sprintf('%s.log',filename),'w+');
        else
            % let's check if we need to update this file
            fid1  = fopen(sprintf('%s.log',filename),'r');
            tline = fgetl(fid1);
            if strcmp(tline,info.date)
                fprintf(1,'Analysis File Up To Date.\n');
                return
            else
                fprintf(1,'Analysis File Changed: Writing Log.\n');
                fclose(fid1);
                fid1 = fopen(sprintf('%s.log',filename),'w+');
            end
        end
    end
else
    % check for pre-existing log file current directory
    loginfo = dir(sprintf('%s.log',filename));
    if isempty(loginfo)
        fprintf(1,'First Use Of Analysis File: Writing Log.\n');
        fid1 = fopen(sprintf('%s.log',filename),'w+');
    else
        % let's check if we need to update this file
        fid1  = fopen(sprintf('%s.log',filename),'r');
        tline = fgetl(fid1);
        if strcmp(tline,info.date)
            fprintf(1,'Analysis File Up To Date.\n');
        else
            fprintf(1,'Analysis File Changed: Writing Log.\n');
            fclose(fid1);
            fid1 = fopen(sprintf('%s.log',filename),'w+');
        end
    end
end
% Open the analysis file to be read
fid2 = fopen(sprintf('%s',[AnalysisFileName '.m']),'r');
fprintf(fid1,'%s\n',info.date);

while ~feof(fid2)
    try
    tline = fgetl(fid2);
    fprintf(fid1,'%s\n',tline);
    catch
        keyboard
    end
end
% Close all files
fclose(fid1);
fclose(fid2);