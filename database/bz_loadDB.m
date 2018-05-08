function [buzsakilab_database] = bz_loadDB(varargin)
% USAGE 
% [buzsakilab_database] = bz_loadDB(varargin)
% 
% INPUTS
%
% DBName       - string, name of database to load
%                [default=buzsakid_metadata]
% 
% forceReload - logical [default=true] When downloading the database the first time, 
%                it will be saved on your local harddrive such that you don't have 
%                download it again. This version will be loaded next time, unless 
%                indicated with the forceReload parameter will indicate that you 
%                want to redownload the database.
% saveMat     - logical [default=false] set to true if you'd like to save
%               DB in a .mat file
% 
% OUTPUT
% buzsakilab_database : a structure with size n x 1, where n is number of sessions. A field for each column 
% in the database.
%
%
% By Peter Petersen
% petersen.peter@gmail.com


p = inputParser;
addParameter(p,'DBname','buzsakid_metadata',@ischar)
addParameter(p,'forceReload',true,@islogical)
addParameter(p,'saveMat',false,@islogical)

parse(p,varargin{:})

DBName = p.Results.DBname;
forceReload = p.Results.forceReload;
saveMat = p.Results.saveMat;

if exist('buzsakilab_database.mat') & forceReload == 0
    disp('Loading local copy of mySQL database')
    load('buzsakilab_database.mat');
else    disp('Loading database from mySQL server')
    % Set preferences
    prefs = setdbprefs('DataReturnFormat');
    setdbprefs('DataReturnFormat','table');
    
    % database connection details
    mySQL_buzsaki.address = '77.104.157.210';
    mySQL_buzsaki.userid = 'buzsakid_users';
    mySQL_buzsaki.password = 'BuzsakiLab9';
    mySQL_buzsaki.database = DBName;
    mySQL_buzsaki.table = 'Datasets';
    % Make connection to database
    conn = database(mySQL_buzsaki.database,mySQL_buzsaki.userid,mySQL_buzsaki.password,'Vendor','MYSQL','Server',mySQL_buzsaki.address,'PortNumber',3306);
    
    % Execute query and fetch results
    curs = exec(conn,['SELECT * FROM ' mySQL_buzsaki.database '.' mySQL_buzsaki.table]);
    if ~strcmp(curs.Message,'Invalid connection.')
        curs = fetch(curs);
        
        buzsakilab_database = table2struct(curs.Data);
        close(curs)
        % Close connection to database
        close(conn)
        % Restore preferences
        setdbprefs('DataReturnFormat',prefs)
        % Clear variables
        clear prefs conn curs
        if saveMat
            save('buzsakilab_database.mat','buzsakilab_database')
        end
        disp('mySQL database loaded successfully')
    else
        disp('Connection failed. mySQL server not available')
        buzsakilab_database = [];
    end
end
