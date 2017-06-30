function version_number = get_ndt_version
% This function returns the neural decoding toolbox version number.
%  It assumes that the current neural decoding toolbox contains a directory
%  with the current toolbox version in the path in the form ndt.x.y.z/


curr_function_dirname = which('get_ndt_version.m');

start_ind = strfind(curr_function_dirname, 'ndt.');


%version_number = strtok(curr_function_dirname(start_ind:end), {'/', '\'});    
   
version_number = strtok(curr_function_dirname(start_ind:end), ['/', '\']);  % changed to be compatible with Octave  


    