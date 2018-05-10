function [bz_datasets,bz_datasets_excluded] = bz_database_submit_collection(database_name,Investigator,NYU_share_path)
% Submit datasets from a folder on the NYU sbare to the Buzsakilab metadata-database
% scans the directory tree for new sessions to add to the database. 
% NYU_Share_Datasets / Investigator / Animal / Session
%
% v1.0
%
% INPUTS
% database_name : 
% A structure with size n x 1, where n is number of sessions to submit.
%
% Investigator : 
% The Investigator level folder to add to the databade
%
% NYU_share_path : 
% Path to the share from your local machine
%
%
% By Peter Petersen
% petersen.peter@gmail.com
cd(fullfile(NYU_share_path,Investigator))
dirAnimals = dir; % dirlist = dirlist(3:end);

bz_datasets = [];
bz_datasets_excluded = [];
k = 1;
kk = 1;
for i = 1:size(dirAnimals,1)
    if ~sum(strcmp(dirAnimals(i).name(1),{'.','_'})) & ~sum(strcmp(dirAnimals(i).name,{'Histology'})) & dirAnimals(i).isdir
        directory_size = folderSizeTree(pwd);
        dirSessions = dir(fullfile(NYU_share_path,Investigator,dirAnimals(i).name));
        for j = 1:size(dirSessions,1)
            if ~sum(strcmp(dirSessions(j).name(1),{'.','_'})) & ~sum(strcmp(dirSessions(j).name,{'Histology'})) & dirAnimals(i).isdir
                
                % Minimum requirements:
                % Folder structure convention: Experimenter / Animal / Session /
                % Files: xml and .dat/.eeg/.lfp and clu,rez,fet
                Session_filebase = fullfile(NYU_share_path,Investigator,dirAnimals(i).name,dirSessions(j).name,dirSessions(j).name);
                if (exist([Session_filebase, '.xml']) & (exist([Session_filebase, '.dat']) | exist([Session_filebase, '.eeg'])| exist([Session_filebase, '.lfp'])))
                    bz_datasets(k).Session = dirSessions(j).name;
                    bz_datasets(k).SessionPath = Session_filebase;
                    bz_datasets(k).Animal = dirAnimals(i).name;
                    bz_datasets(k).Investigator = Investigator;
                    bz_datasets(k).Sessions = 1;
                    bz_datasets(k).UsageLicence = 'Publicly available';
                    bz_datasets(k).NYUSharefoldersizes = directory_size.size{find(strcmp(fullfile(NYU_share_path,Investigator,dirAnimals(i).name,dirSessions(j).name),directory_size.name))}/(1024*1024);
                    bz_datasets(k).NYUWebsharePath = ['https://buzsakilab.nyumc.org/datasets/', bz_datasets(k).Investigator,'/',bz_datasets(k).Animal,'/',bz_datasets(k).Session,'/||Download ' bz_datasets(k).Session];
                    bz_datasets = bz_database_extract_meta(bz_datasets,k);
                    fprintf(['Dataset size: ', num2str(bz_datasets(k).NYUSharefoldersizes,'%10.0f'),'MB. '])
                    k = k + 1;
                else
                    bz_datasets_excluded(kk).Session = dirSessions(j).name;                    
                    bz_datasets_excluded(kk).SessionPath = Session_filebase;
                    bz_datasets_excluded(kk).Animal = dirAnimals(i).name;
                    bz_datasets_excluded(kk).Investigator = Investigator;
                    kk = kk + 1;
                end
            end
        end
    end
end

fprintf(['\n',num2str(size(bz_datasets,2)), ' sessions detected. ', num2str(size(bz_datasets_excluded,2)), ' sessions not fulfilling requirements'])
bz_datasets = rmfield(bz_datasets,{'SessionPath'});
bz_database_submit(bz_datasets,database_name);