function buzsakilab_database = load_buzsakilab_database_datasets(force_reload)
% Importing Data from the Buzsakilab metadata-database
% v1.1
% 
% INPUT
% force_reload : When downloading the database the first time, 
% it will be saved on your local harddrive such that you don't have 
% download it again. This version will be loaded next time, unless 
% indicated with the force_reload parameter will indicate that you 
% want to redownload the database.
% 
% OUTPUT
% buzsakilab_database : a structure with a field for each column 
% in the database.
%
%
% By Peter Petersen
% petersen.peter@gmail.com

if nargin == 0
    force_reload = 0;
end

if exist('buzsakilab_database.mat') & force_reload == 0
    disp('Loading local copy of mySQL database')
    load('buzsakilab_database.mat');
else
    disp('Loading database from mySQL server')
    % Set preferences
    prefs = setdbprefs('DataReturnFormat');
    setdbprefs('DataReturnFormat','structure');
    
    % Make connection to database
    mySQL_buzsaki.database = 'buzsakid_metadata';
    mySQL_buzsaki.userid = 'buzsakid_contrib';
    mySQL_buzsaki.password = 'BuzsakiLab2017';
    mySQL_buzsaki.table = 'Datasets';
    conn = database(mySQL_buzsaki.database,mySQL_buzsaki.userid,mySQL_buzsaki.password,'Vendor','MYSQL','Server','77.104.157.210','PortNumber',3306);
    
    % Execute query and fetch results
    curs = exec(conn,['SELECT * FROM ' mySQL_buzsaki.database '.' mySQL_buzsaki.table]);
    if ~strcmp(curs.Message,'Invalid connection.')
        curs = fetch(curs);
        buzsakilab_database = curs.Data;
        close(curs)
        % Close connection to database
        close(conn)
        % Restore preferences
        setdbprefs('DataReturnFormat',prefs)
        % Clear variables
        clear prefs conn curs
        save('buzsakilab_database.mat','buzsakilab_database')
        disp('mySQL databse loaded successfully')
    else
        disp('Connection failed. mySQL server not available')
        buzsakilab_database = [];
    end
end
