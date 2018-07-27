Welcome to the Buzsaki lab repository! The goal is to have this repo as your one-stop-shop for code you may need as a Buzsaki Lab member. This will include preprocessing pipelines as well as analysis functions. 
IMPORTANT: everything is under collective development by the lab. If you find something that needs changing/adding, please open an issue and improve it (:


New here?
FIRST: Run the function 'compileBuzcode' in the top directory.
SECOND: Check out the tutorials in the 'tutorials' folder.
THIRD: Check out the wiki. This includes data formatting standards, instructions on using git, and other stuff. As you find things that aren't clear, consider improving the wiki to make them more clear for the next user.




DATA FORMATTING STANDARDS
Everything will conform to a single (flexible, but documented) database structure. All files pertaining to a single recording will live in a single self-contained folder, or a 'basePath' (whateverPath/baseName/), where baseName is the name of the recording. Files will follow the naming convention basePath/baseName.fileName.filetype. For example, /recording7/recording7.ripples.events.mat will be a file containing information about ripples from a recording named recording7, it’s contents will be in the format prescribed to .events.mat type files. 

More on this can be found in the wiki here:
https://github.com/buzsakilab/buzcode/wiki/Data-Formatting-Standards

As you (inevitably) run into data types that don't really fit into any of these boxes, please open an issue to discuss with the group, and add the necessary format/standards to the wiki.



