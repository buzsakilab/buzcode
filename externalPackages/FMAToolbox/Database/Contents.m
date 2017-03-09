% Database and batch processing functions for FMAToolbox.
%
% These functions can be used to perform batch data processing and store data
% and figures into one or more databases.
%
% Batch Processing
%
% Batch processing is useful if you need to repeatedly perform a given analysis
% on different data sets. This can be achieved using <a href="StartBatch">StartBatch</a>, which will
% repeatedly run a given 'batch' function on each data set. <a href="StartBatch">StartBatch</a> provides
% a number of advanced features such as error resilience, error logging, delayed
% processing, etc. The 'batch' function typically performs computations and
% stores the results in a database using <a href="DBAddVariable">DBAddVariable</a> and/or <a href="DBAddFigure">DBAddFigure</a>. Once
% the batch is completed, the results can be retrieved from the database using
% <a href="DBGetValues">DBGetValues</a>. Different sets of results can be combined using <a href="DBMatchValues">DBMatchValues</a>.
%
% Database Management
%
% Databases can be private and run on the local machine, or central databases
% running on a remote server and shared between many different users. Data and
% figures can be tagged (experiment ID, name, comments) so you can later search
% specific data subsets. Programs and parameters used to generate the data can
% also be stored for reference. Figures can be stored as matlab fig files
% and/or PNG images, and be later exported as HTML galleries which are easy to
% share with other investigators.
%
% Notes
%
%  The database backend is MySQL and the interface is provided by the <a href="http://sourceforge.net/projects/mym/">mYm</a>
%  toolbox by Y. Maret (EPFL, Lausanne, Switzerland).
%
%  Although not required, installing <a href="http://sites.google.com/site/oliverwoodford/software/export_fig">export_fig</a> by O. Woodford (Oxford
%  University, England) will yield better PNG images.
%
% Batch processing
%
%   StartBatch           - Start a new batch job.
%   ShowBatch            - Show data sets in a batch job.
%   GetBatch             - Get batch job output.
%   BatchInfo            - Get batch job information.
%   CancelBatch          - Cancel batch job.
%   CleanBatches         - Delete completed batch jobs from memory.
%   DebugBatch           - Assign variables to help debug a batch job.
%   Store                - Store variable in memory to avoid redundant computations.
%   Recall               - Recall variable from memory to avoid redundant computations.
%
% Database management.
%
%   DBCreate             - Create a new database (requires admin rights)
%   DBCreateTables       - Create tables for variables and figures
%   DBConnect            - Connect to the database server
%   DBUse                - Set (or determine) current database
%   DBList               - List existing databases
%   DBDuplicate          - Duplicate database.
%   DBMerge              - Merge databases.
%   DBRemove             - Remove database.
%
% Storing data/figures.
%
%   DBAddFigure          - Store a figure (+ tags, programs, etc.)
%   DBAddVariable        - Store a variable (+ tags, programs, etc.)
%   DBRemoveFigures      - Remove all figures that match given criteria.
%   DBRemoveVariables    - Remove all variables that match given criteria.
%
% Searching data/figures.
%
%   DBGetFigures         - Get all figures that match given criteria
%   DBGetValues          - Group values for all variables that match given criteria
%   DBMatchValues        - Match values for two different sets of variables
%   DBGetVariables       - Get all variables that match given criteria
%   DBListFigures        - List figures matching certain criteria.
%   DBListVariables      - List variables matching certain criteria.
%
% Exporting figures.
%
%   DBExportGallery      - Create an HTML gallery from a set of figures
%