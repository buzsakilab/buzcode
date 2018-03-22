function [newDB] = bz_updateSQLFromFiles(nyuSharePath)
% USAGE
% [newDB] = bz_updateSQLFromFiles()
%
% INPUTS
%
%
%
%
%
% INPUTS
%
%
%
% DESCRIPTION
% This function takes a directory path as input (from /Datasets/ folder) 
% and searches for recording session folders to add to the SQL database
%
% david tingley, 2018


experimenterList = dir([nyuSharePath '/Buzsakilabspace/Datasets/*']);

totalSessions = 0;
%% step 1, find all sessions and compile a path list
for exp = 1:length(experimenterList)
    warning off
    if ~strcmp(experimenterList(exp).name,'unsorted') & ~strcmp(experimenterList(exp).name(1),'.') % ignore unsorted folder and hidden files
        
        experimenterList(exp).name
        [isBuzcode{exp}, buzcodePaths{exp}] = dirwalk([nyuSharePath '/Buzsakilabspace/Datasets/' experimenterList(exp).name],@bz_isBuzcode,3);
        % finds recordings that are currently in buzcode format
        
        [isSession{exp}, sessionPaths{exp}] = dirwalk([nyuSharePath '/Buzsakilabspace/Datasets/' experimenterList(exp).name],@bz_isSession,3);
        % finds recordings that have a minimal dataset (*.dat or *.eeg file)
        
        bzSessions(exp) = sum(cell2mat(isBuzcode{exp}));
        totalSessions(exp) = sum(cell2mat(isSession{exp}));
        totalSessions = totalSessions + sum(cell2mat(isSession))
    end
end

%% step 2, get current DB and check if these already exist 

    db = bz_loadDB('dbname','buzsakid_submitdataset');
    
%% step 3, add new recordings

