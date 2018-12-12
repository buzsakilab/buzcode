function bz_database = bz_database_load(database_name,filters)
% Loading the Buzsakilab metadata-database
%
% v1.3
%
% INPUTS
% load_database : 
% Which database to load. Default: 'buzsakid_metadata'
%
% filters : 
% Filters by column names can be applied according to mySQL standards:
% Conditions: AND, OR
% Operators: =, !=, >, <, >=, <=, LIKE
% Parentheses () AND ()
% List of values: IN 
%
% Eg. 'Animal = "Cicero"'
%     'Species != "Rat"'
%     'Id  1' 
%     'BrainRegion IN ("Hippocampus","Thalamus")'
%     'Sex = "Male" AND Species = "Rat"'
%
%
% OUTPUT
% bz_database : 
% A structure with size n x 1, where n is number of sessions. A field for each column
% in the database with optionally applied filterss.
%
%
% By Peter Petersen
% petersen.peter@gmail.com

if ~exist('database_name','var') || isempty(database_name)
    database_name = 'buzsakid_metadata';
end
if ~exist('filters','var')
    filters = 0;
end

fprintf(['\nLoading database: ', database_name, ' from server. '])

% Set preferences
prefs = setdbprefs('DataReturnFormat');
setdbprefs('DataReturnFormat','table');

% Make connection to database
bz_mySQL = bz_database_credentials; % Loading database credentials

conn = database(database_name,bz_mySQL.userid,bz_mySQL.password,'Vendor','MYSQL','Server',bz_mySQL.address,'PortNumber',3306);

% Execute query and fetch results
if ~ischar(filters)
    curs = exec(conn,['SELECT * FROM ' database_name '.' bz_mySQL.table]);
else
    disp(['Filter applied: ' filters])
    curs = exec(conn,['SELECT * FROM ' database_name '.' bz_mySQL.table ' WHERE ' filters]);
end
if ~strcmp(curs.Message,'Invalid connection.')
    curs = fetch(curs);
    
    bz_database = table2struct(curs.Data);
    close(curs)
    close(conn) % Close connection to database
    setdbprefs('DataReturnFormat',prefs) % Restore preferences
    clear prefs conn curs % Clear variables
    disp([num2str(size(bz_database,1)),' datasets loaded successfully from the database'])
else
    disp('Connection failed. mySQL server not available.')
    bz_database = [];
    warning(curs.Message)
end
