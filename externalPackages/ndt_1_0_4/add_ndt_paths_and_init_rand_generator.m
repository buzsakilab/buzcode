function add_ndt_paths_and_init_rand_generator

% This function adds all the appropriate directories for the NDT toolbox to work
%   and initializes the random number generator based on the clock (by default
%   Matlab uses the same seed everything it is started which can lead to unwanted behavior).

%==========================================================================

%     This code is part of the Neural Decoding Toolbox.
%     Copyright (C) 2011 by Ethan Meyers (emeyers@mit.edu)
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
%==========================================================================   





% file the path where the toolbox is 
toolbox_basedir_name = [fileparts(which('add_ndt_paths_and_init_rand_generator')) '/']

% add all the appropriate sub-directories
addpath([toolbox_basedir_name 'datasources/']);
addpath([toolbox_basedir_name 'feature_preprocessors/']);
addpath([toolbox_basedir_name 'classifiers/']);
addpath([toolbox_basedir_name 'cross_validators/']);
addpath([toolbox_basedir_name 'tools/']);
addpath([toolbox_basedir_name 'helper_functions/']);
addpath(genpath([toolbox_basedir_name 'external_libraries/']));  % include subdirectories of any external libraries   

if isOctave
    addpath([toolbox_basedir_name 'octave_code/']);
    'Running Octave - adding Octave specific code to the search path'
end
    
% initialize the random number generator
rand('state',sum(100*clock)); 
disp('initializing the matlab random number generator to an aribrary clock value, i.e., the ''rand'' function should now work properly');
















