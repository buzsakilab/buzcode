% Reading mySQL database
if ispc
    NYU_share_path =  'Q:\Buzsakilabspace\Datasets';
elseif ismac
    NYU_share_path = '/Volumes/buzsakilab/Buzsakilabspace/Datasets';
end
database_name = 'buzsakid_metadata'; % buzsakid_metadata, buzsakid_submitdataset
bz_database = bz_database_load(database_name,'BrainRegion LIKE "%Hippocampus%"');

%% % Submitting dataset to database
bz_datasets = bz_database([18:90]);
bz_database_submit(bz_datasets,'buzsakid_submitdataset');

%% % Example script for loading and plotting a dataset from the database
j = 1;
session_path = fullfile(NYU_share_path, bz_database(j).Investigator, bz_database(j).Animal, bz_database(j).Session);
cd(session_path)
xml = LoadXml(fullfile(session_path,[bz_database(j).Session,'.xml']));
units = loadClusteringData(bz_database(j).Session,bz_database(j).Unitsformat,session_path,1);
% Plotting a histogram of the firing rates across sorted units
figure, hist(horzcat(units.total)), title('Histogram of spikes/unit')

%% % Pulling out information from the dataset folders
bz_database = bz_database_load;
datasets = [1:size(bz_database,1)];
bz_datasets = bz_database;
for j = datasets
    bz_datasets(j).SessionPath = fullfile(NYU_share_path,bz_datasets(j).Investigator,bz_datasets(j).Animal, bz_datasets(j).Session,bz_datasets(j).Session);
    bz_datasets = bz_database_extract_meta(bz_datasets,j);
end
bz_datasets = rmfield(bz_datasets,{'SessionPath'});
bz_database_update(bz_datasets(datasets),'buzsakid_metadata');% buzsakid_metadata,buzsakid_submitdataset

%% % Updating an existing dataset in database
bz_database_update(bz_database(1));

%% % Scanning a directory for new sessions
Dataset_directory = 'BerenyiT'; % BerenyiT
database_name = 'buzsakid_metadata'; % buzsakid_metadata,buzsakid_submitdataset
[bz_datasets,bz_datasets_excluded] = bz_database_submit_collection(database_name,Dataset_directory,NYU_share_path);
