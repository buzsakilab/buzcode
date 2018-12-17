function buzsakilab_database_submit(bz_datasets,database_name)
% Update existing dataset in the Buzsakilab metadata-database
%
% v1.3
%
% INPUTS
% bz_datasets :
% structure with subfields named after the mySQL columns. The subfields has
% to have the right format to make updates to the database.
%
% database_name :
% Name of database to submit to. Default 'buzsakid_submitdataset'
%
% By Peter Petersen
% petersen.peter@gmail.com

if nargin == 1
    database_name = 'buzsakid_submitdataset';
end

fprintf(['\nUpdating dataset(s) in database: ', database_name, ' from server. \n'])

% Make connection to database
bz_mySQL = bz_database_credentials; % loading credentials

conn = database(database_name,bz_mySQL.userid,bz_mySQL.password,'Vendor','MYSQL','Server',bz_mySQL.address,'PortNumber',3306);
setdbprefs('NullNumberWrite','NaN');
setdbprefs('DataReturnFormat','table');
curs = exec(conn,['SELECT * FROM ' database_name '.' bz_mySQL.table]);
curs = fetch(curs);

if ~strcmp(curs.Message,'Invalid connection.')
    bz_database = table2struct(curs.Data);
    % attributes = attr(curs);
    close(curs);
    
    for j = 1:length(bz_datasets)
        bz_dataset = bz_datasets(j);
        id = find(bz_dataset.Id==[bz_database.Id]);
        if ~isempty(id) && strcmp(bz_dataset.Session,bz_database(id).Session) & strcmp(bz_dataset.Animal,bz_database(id).Animal) & strcmp(bz_dataset.Investigator,bz_database(id).Investigator)
            
            whereclause = ['WHERE Id = ' num2str(bz_dataset.Id)];
            disp(['Updating dataset ' num2str(bz_dataset.Id) ': ' bz_dataset.Investigator, ', ', bz_dataset.Animal,', ', bz_dataset.Session])
            bz_dataset = rmfield(bz_dataset, 'Id');
            
            if any(structfun(@isempty, bz_dataset))
                idx = find(structfun(@isempty, bz_dataset));
                names = fieldnames(bz_dataset);
                for field = 1:length(idx)
                    % warning([names{idx(field)} ' was empty, filling with null.'])
                    bz_dataset = setfield(bz_dataset,names{idx(field)},'null');
                end
            end
            
            % Removing fields with NaN and null values
            data_table = struct2table(bz_dataset);
            temp = fieldnames(bz_dataset);
            temp2 = find(ismissing(data_table,{NaN,'null'}));
            bz_dataset = rmfield(bz_dataset,temp(temp2));
            % Data to insert
            colnames = fieldnames(bz_dataset)'; %
            data_table = struct2table(bz_dataset);
            
            % data_table = buzsakilab_database_verify_fieldtypes_and_dataformat(data_table,attributes)
            update(conn,bz_mySQL.table,colnames,data_table, whereclause);
        else
            warning(['Dataset with id=' num2str(bz_dataset.Id) ' not consistent with Investigator, Animal and Sessions in database with same id'])
        end
    end
else
    disp('Connection failed. mySQL server not available')
end
% Close connection to database
close(conn)
% Clear variables
clear conn curs
