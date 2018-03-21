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
for exp = 1:length(experimenterList)
    if ~strcmp(experimenterList{exp},'unsorted') % ignore unsorted folder
        
        dirwalk([nyuSharePath '/Buzsakilabspace/Datasets/' experimenterList{exp}],@bz_isBuzcode)
        % finds recordings that are currently in buzcode format
        
%         dirwalk([nyuSharePath '/Buzsakilabspace/Datasets/' experimenterList{exp}],@bz_isSession)
        % finds recordings that have a minimal dataset (*.dat or *.eeg
        % file)
        
        
        
    end
end