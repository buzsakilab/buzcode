function buzsakilab_database_submit(bz_datasets,database_name)
% Submit dataset to the Buzsakilab metadata-database
%
% v1.1
%
% INPUT
% bz_datasets : 
% A structure with size n x 1, where n is number of sessions to submit.
% 
% database_name :
% Name of database to submit to. Default 'buzsakid_submitdataset'
%
%
% By Peter Petersen
% petersen.peter@gmail.com

if nargin == 1
    database_name = 'buzsakid_submitdataset';
end

fprintf(['\nSubmitting datasets to database: ''', database_name, '''\n'])

% Make connection to database
bz_mySQL = bz_database_credentials; % Credentials loaded

conn = database(database_name,bz_mySQL.userid,bz_mySQL.password,'Vendor','MYSQL','Server',bz_mySQL.address,'PortNumber',3306);
curs = exec(conn,['SELECT * FROM ' database_name '.' bz_mySQL.table]);
curs = fetch(curs);

if ~strcmp(curs.Message,'Invalid connection.')
    bz_database = table2struct(curs.Data);
    % attributes = attr(curs);
    close(curs);
    
    % Execute query and fetch results
    if isfield(bz_datasets, 'Id')
        bz_datasets = rmfield(bz_datasets, 'Id');
    end
    
    for k = 1:size(bz_datasets,1)
        if sum(strcmp(bz_datasets(k).Session,{bz_database.Session}) & strcmp(bz_datasets(k).Animal,{bz_database.Animal}) & strcmp(bz_datasets(k).Investigator,{bz_database.Investigator}))
            bz_datasets(k).NewDataset = 0;
        else
            bz_datasets(k).NewDataset = 1;
        end
    end
    clear bz_database;
    n_datasets = size(bz_datasets,1);
    bz_datasets = bz_datasets(find([bz_datasets.NewDataset]))';
    
    if size(bz_datasets,2) > 0
        bz_datasets = rmfield(bz_datasets,{'NewDataset'});
        for j = 1:length(bz_datasets)
            % Checking for existing datasets in the database
            disp([num2str(j) '. dataset: ' bz_datasets(j).Investigator, ', ', bz_datasets(j).Animal,', ', bz_datasets(j).Session])
            if any(structfun(@isempty, bz_datasets(j)))
                idx = find(structfun(@isempty, bz_datasets(j)));
                names = fieldnames(bz_datasets);
                for field = 1:length(idx)
                    warning([names{idx(field)} ' was empty, filling with null.'])
                    bz_datasets(j) = setfield(bz_datasets(j),names{idx(field)},'null');
                end
            end
            % Data to insert
            colnames = fieldnames(bz_datasets(j))'; %
            data_table = struct2table(bz_datasets(j));
            % data_table = buzsakilab_database_verify_fieldtypes_and_dataformat(data_table,attributes)
            
            insert(conn,bz_mySQL.table,colnames,data_table);
        end
        disp([num2str(size(bz_datasets,2)), '/' num2str(n_datasets), ' dataset(s) submitted to database. ', num2str(n_datasets-size(bz_datasets,2)) ' dataset(s) already exist in database'])        
    else
        disp('Datasets already exist in database')
    end
else
    disp('Connection failed. mySQL server not available')
    disp(curs.Message)
end
close(conn) % Close connection to database
clear conn curs % Clear variables
