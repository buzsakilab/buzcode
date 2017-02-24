% Database access functions for FMAToolbox.
%
% These functions can be used to store data and figures into a database. The
% database can be a private database running on the local machine, as well
% as a central database running on a remote server and shared between many
% different users. Data and figures can be tagged (experiment ID, item name,
% comments). Programs and parameters used to generate them can also be stored
% for reference. Database queries (in SQL) can then be used to search for
% specific sets of data.
%
% Figures are stored as matlab fig files and/or PNG images. Figures can be
% exported to HTML galleries to easily share them with other investigators.
%
% Notes
%
%  The database backend is MySQL and the interface is provided by the <a href="http://sourceforge.net/projects/mym/">mYm</a>
%  toolbox by Y. Maret (EPFL, Lausanne, Switzerland).
%
%  Although not required, installing <a href="http://sites.google.com/site/oliverwoodford/software/export_fig">export_fig</a> by O. Woodford (Oxford
%  University, England) will yield better PNG images.
%
% Database management.
%
%   DBCreate             - Create a new database (requires admin rights)
%   DBCreateTables       - Create tables for variables and figures
%   DBConnect            - Connect to the database server
%   DBUse                - Use an existing database
%
%
% Storing data/figures.
%
%   DBInsertFigure       - Store a figure (+ tags, programs, etc.)
%   DBInsertVariable     - Store a variable (+ tags, programs, etc.)
%
% Searching data/figures.
%
%   DBSelectFigure       - Find a uniquely identified figure
%   DBSelectFigures      - Find all figures that match given criteria
%   DBSelectVariable     - Find a uniquely identified variable
%   DBSelectVariables    - Find all variables that match given criteria
%
% Listing data/figures.
%
%   DBListFigures        - List all figures in the database
%   DBListVariables      - List all variables in the database
%
% Exporting figures.
%
%   DBExportGallery      - Create an HTML gallery from a set of figures
%