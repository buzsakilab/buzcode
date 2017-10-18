Welcome to the Buzsaki lab repository!
The goal is to have this repo as your one-stop-shop for code you may need as a Buzsaki Lab member. This will include preprocessing pipelines as well as analysis functions. Everything will conform to a single (flexible, but documented) database structure.

IMPORTANT: everything is under collective development by the lab.

DATA FORMATTING STANDARDS
(under development: read and contribute to the wiki!)
https://github.com/buzsakilab/buzcode/wiki/Data-Formatting-Standards

The basic philosophy is built on Neuroscope data formatting, with supplementary .mat types for commonly used data types, each which have a standardized format as outlined in the git wiki. All files pertaining to a given recording should be in a single self-contained folder called baseName. 
For example, /recording7/recording7.ripples.events.mat will be a file containing information about ripples from a recording named recording7, it’s contents will be in the format prescribed to .events.mat type files. 

As you (inevitably) run into data types that don't really fit into any of these boxes, please consult the #code-development channel at buzsakilab.slack.com/messages/code-development/, and add the necessary format/standards to the wiki.






