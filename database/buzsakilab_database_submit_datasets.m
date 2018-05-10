function buzsakilab_database_submit_datasets(buzsakilab_datasets)
% Submit dataset to the Buzsakilab metadata-database
% v1.0
%
% INPUT
% buzsakilab_datasets
% structure with subfields named after the mySQL columns.
%
%
% By Peter Petersen
% petersen.peter@gmail.com
if length(buzsakilab_datasets) > 1
    disp('Submitting datasets to the buzsakilab database')
else
    disp('Submitting dataset to the buzsakilab database')
end

% Getting attributed for the fields
mySQL_buzsaki.address = '77.104.157.210';
mySQL_buzsaki.userid = 'buzsakid_users';
mySQL_buzsaki.password = 'BuzsakiLab9';
mySQL_buzsaki.database = 'buzsakid_submitdataset';
mySQL_buzsaki.table = 'Datasets';

% Make connection to database
conn = database(mySQL_buzsaki.database,mySQL_buzsaki.userid,mySQL_buzsaki.password,'Vendor','MYSQL','Server',mySQL_buzsaki.address,'PortNumber',3306);
curs = exec(conn,['SELECT * FROM ' mySQL_buzsaki.database '.' mySQL_buzsaki.table]);
curs = fetch(curs);

if ~strcmp(curs.Message,'Invalid connection.')
    attributes = attr(curs);
    close(curs);
    close(conn);
    
    % database connection details
    mySQL_buzsaki.database = 'buzsakid_submitdataset';
    conn = database(mySQL_buzsaki.database,mySQL_buzsaki.userid,mySQL_buzsaki.password,'Vendor','MYSQL','Server',mySQL_buzsaki.address,'PortNumber',3306);
    
    % Execute query and fetch results
    if isfield(buzsakilab_datasets, 'Id')
        buzsakilab_datasets = rmfield(buzsakilab_datasets, 'Id');
    end
    for j = 1:length(buzsakilab_datasets)
        disp(['Submitting dataset ' num2str(j) ': ' buzsakilab_datasets(j).Session])
        %% check for empty fields, and will with 'null'
        if any(structfun(@isempty, buzsakilab_datasets(j)))
            idx = find(structfun(@isempty, buzsakilab_datasets(j)));
            names = fieldnames(buzsakilab_datasets);
            for field = 1:length(idx)
                warning([names{idx(field)} ' was empty, filling with null.'])
                buzsakilab_datasets(j) = setfield(buzsakilab_datasets(j),names{idx(field)},'null');
            end
        end
        %% check if this session already exists in the DB
        
        
        %% Data to insert
        colnames = fieldnames(buzsakilab_datasets(j))'; %
        data_table = struct2table(buzsakilab_datasets(j));
        % data_table = buzsakilab_database_verify_fieldtypes_and_dataformat(data_table,attributes)
        insert(conn,mySQL_buzsaki.table,colnames,data_table);
    end
    
    % Close connection to database
    close(conn)

    % Clear variables
    clear conn curs
    if length(buzsakilab_datasets) > 1
        disp('Datasets submitted successfully')
    else
        disp('Dataset submitted successfully')
    end
else
    disp('Connection failed. mySQL server not available')
end
