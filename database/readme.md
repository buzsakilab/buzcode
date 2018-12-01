# Meta-data database
__The database is built in mySQL and contains metadata describing datasets collected by Investigators from the Buzsaki lab hosted on the share NYU. The following functions can be used to interact with the database:__

* bz_datasbase_load: Loading the database into Matlab. Filters can be applied, using mySQL standards (see the documentation of the file for examples)<br>
* bz_database_update: Update an existing dataset in the database<br>
* bz_datasbase_submit: Submitting a dataset  to the database<br>
* bz_database_submit_collection: Submitting a dataset collection to the database provide investigator_name, database_name, path the share as input.<br>
* bz_database_credentials: the credentials used when interacting with the database. The provided user account has reading rights to the main database (buzsakid_metadata), yet datasets can still be submitted to the secondary database (buzsakid_submitdataset) <br>
* bz_database_example_scripts: a list of example calls to the database using above functions.

The NYU share dataset folder is organized across three levels: Investigator/Animal/Session. Folders beginning with . or _ and reserverd names are not scanned for electrophysiological data (at this point only Histology).

### Webshare
You can access the datasets from our webshare: [buzsakilab.nyumc.org/datasets/](https://buzsakilab.nyumc.org/datasets/).

### Web interface
The database also has a webinterface that can be used on the [Buzsakilab website](http://buzsakilab.com/wp/datasets/).

### Matlab installation instructions and more
Learm more about [installing the mySQL database tool and setup Matlab here](http://buzsakilab.com/wp/database-access/).
