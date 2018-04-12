function mypath = install
%INSTALL Set search paths, and download and install dependencies.
% 
%   IOSR.INSTALL downloads and installs the toolbox dependencies. The
%   function also adds the required paths to the Matlab search path.
% 
%   MYPATH = IOSR.INSTALL returns the old Matlab search path MYPATH.

%   Copyright 2016 University of Surrey.
    
    %% download and install SOFA

    % install dir
    currdir = cd;
    cd([fileparts(which(mfilename('fullpath'))) filesep '..']);
    directory = pwd;
    sofa_folder = [directory filesep 'deps' filesep 'SOFA_API'];
    
    if ~(exist(sofa_folder,'dir') == 7)
        % download and install
        sofa_filename = 'sofa-api.zip';
        try % Sourceforge location changes from time to time
            websave(sofa_filename,'http://vorboss.dl.sourceforge.net/project/sofacoustics/sofa-api-mo-1.0.2.zip');
        catch % Fall back to development release on GitHub
            display('Warning: Failed to download SOFA v1.0.2 from Sourceforge, downloading development release...')
            websave(sofa_filename,'https://github.com/sofacoustics/API_MO/archive/master.zip');
        end
        unzip(sofa_filename,sofa_folder);
        movefile([sofa_folder filesep 'API_MO' filesep '*'],[sofa_folder filesep]);
        
        % clean up
        delete(sofa_filename)
        rmdir([sofa_folder filesep 'doc' filesep],'s')
        rmdir([sofa_folder filesep 'HRTFs' filesep],'s')
        rmdir([sofa_folder filesep 'API_MO' filesep],'s')
    else
        display(strcat('Found existing SOFA directory: ', sofa_folder)) 
    end
    
    %% Add directories to path
    
    cd(directory);
    mypath = addpath(directory,...
        [directory filesep 'deps' filesep 'SOFA_API']);
    
    %% start SOFA
    
    SOFAstart(0);
    
    cd(currdir); % return to original directory
    
end
