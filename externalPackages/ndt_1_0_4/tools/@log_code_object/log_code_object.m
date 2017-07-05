classdef log_code_object < handle

% This helper object saves the code of files that have been run. It is useful for saving information about
%  the exact code that was run so that one has a record of how particular results were generated.  This object
%  gives several different methods for saving code in different files/memory configurations.  The results  are in a structure that
%  has the name and path of the file being recorded, the text in the body of the file, the time when the file was last
%  modified and the time when the log_code_object recorded the file This method inherits from the from the handle
%  class so all logged files are persistently kept in the logged_code_structure as new files are added.   
%
%  The constructor for this object takes no parameters, i.e., to create an object use:  
%       log_code_obj = log_code_object; 
%
%  The methods that can be used to log particular files are:
%
%   1.  log_current_file(log_code_obj).  This method logs the file that the called this method.  
%
%   2.  log_specfic_files(log_code_obj, file_names).  If file_name is a string, this method logs
%         the file that is in file_names, or if files_names is a cell array of strings, this method
%         logs all the files that are in the cell array.  
%   
%   3.  log_files_in_directory(log_code_obj, directory_name, file_pattern).  This method logs all the files
%         the are in the directory directory_name that have the pattern listed in file_pattern.  For example,
%         to log all the m-files in the current directory one would use the methods:  
%         log_code_obj.log_files_in_directory('./', '*.m')  
%
%   4. log_functions_in_memory(log_code_obj, directory_name).  This method logs all the functions that are in 
%       memory from a particular directory given in directory_name.  This is particular useful if there are a
%       large number of files in a given directory but you only want to save the ones that were used by your code
%       (see examples below).
%
%
%  To get the logged information one can call the following method (or just access the field logged_code_structure):
%
%   logged_code_structure = return_logged_code_structure(log_code_obj)  
%
%   The number of elements in logged_code_structure is equal to the number of files logged, and the structure has the following fields:  
%
%       1. .code_filename:  the name of the file that has been recorded 
%       2. .code_filepath:  the full path to the directory where the file has been recorded 
%       3. .code_body:  the text that is in the file
%       4. .last_modified_time:  the time when the file was last modified
%       5. .logged_time:  the time the file was logged (useful in case the file was changed after being logged).
%
%
%   Examples of use:
%
%    We recommend calling calling the method log_code_obj.log_current_file   in the beginning of
%       any code you have the specifies what decoding object/parameters are used in a data analysis (that way
%       if you change the code while the decoding analysis is running the original parameters will be recorded).  
%       If one wants to log which functions in the Neural Decoding Toolbox were used, one can then call the 
%       log_code_obj.log_functions_in_memory(path_to_toolbox) method to save all the functions that were used.
%         
%       For example:
%
%           % start of decoding scripts
%
%           log_code_obj = log_code_object; 
%           log_code_obj.log_current_file;   % log the script that is currently being run
%           
%           % create the decoding objects, run the decoding analysis with DECODING_RESULTS = cv.run_cv_decoding;
%   
%           log_code_obj.log_functions_in_memory(path_to_toolbox);  % log the functions that were used in the NDT
%
%           decoding_code_log = log_code_obj.logged_code_structure;
%           save('my_results', 'DECODING_RESULTS', 'decoding_code_log');
%       
%          

% The code was modified from some functions written by Jim Mutch.

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



    
    
    properties 
    
        % the main structure that contains all the information recorded by this object
        logged_code_structure = struct('code_filename', {}, 'code_filepath', {}, 'code_body', {}, 'last_modified_time', {}, 'logged_time', {});

    end
    
    
    
    methods
        
        function log_current_file(log_code_obj)
        % This function adds the code that called this method to the logged_code_structure
            
            dbstack_structure = dbstack('-completenames');
                        
            if numel(dbstack_structure) == 1
                warning('one needs to call the method log_current_file from a .m file (calling it from the command line does not make sense) so no files have been logged')
            else
                current_caller_full_file_name = dbstack_structure(2).file;           
                log_code_obj.add_file_to_logged_code_structure(current_caller_full_file_name);
            end         
                
        end
        
        
        
        function log_functions_in_memory(log_code_obj, directory_name)
        %  This function adds all the functions in memory that are in the directory named 'directory_name' to the logged_code_structure  
            
            all_function_names_in_memory = inmem('-completenames');

            function_names_in_directory_that_are_in_memory = all_function_names_in_memory(strmatch(directory_name, all_function_names_in_memory));

            if isempty(function_names_in_directory_that_are_in_memory), warning('There are no function names in directory_name that are currently in memory, so no files have been added to the logged_code_structure'); end
                
            for iFile = 1:numel(function_names_in_directory_that_are_in_memory)
                log_code_obj.add_file_to_logged_code_structure(function_names_in_directory_that_are_in_memory{iFile});
            end
                   
        end
        
        
        
         function log_files_in_directory(log_code_obj, directory_name, file_pattern)
        %  This function adds all the function in the directory named 'directory_name' to the logged_code_structure   
               
            the_dir = dir([directory_name file_pattern]);
        
            if isempty(the_dir), warning('There are no files with the specific file_pattern in the current directory, so no files have been added to the logged_code_structure'); end
            
            for iFile = 1:numel(the_dir)
                log_code_obj.add_file_to_logged_code_structure(the_dir(iFile).name);
            end

         end
        
        
        function log_specfic_files(log_code_obj, file_names)
        % function_names can be a string in which the function given in the string name will be logged, 
        %  or it can be a cell array of strings in which all functions names given in the cell array will be logged.
           
            if isstr(file_names)
            
                log_code_obj.add_file_to_logged_code_structure(file_names);
            
            elseif iscell(file_names)
                
                for iFile = 1:numel(file_names)
                    log_code_obj.add_file_to_logged_code_structure(file_names{iFile});
                end
                
            else
                 warning('file_names must be a string with a file name to be added to the logged_code_structure, or a cell array containing the names that should be added, No files have been added to the logged_code_structure'); 
            end
            
        end


        function logged_code_structure = return_logged_code_structure(log_code_obj)
        % This function returns the logged code structure
            logged_code_structure = log_code_obj.logged_code_structure;            
        end
        
        
    end
    
    
    
    
    methods (Access = 'private')
        
 
        function add_file_to_logged_code_structure(log_code_obj, full_file_name)
        % Helper function that adds information to the logged_code_structure
            
            % index in the logged_code_structure where the file should be added
            curr_log_code_struct_ind = numel(log_code_obj.logged_code_structure) + 1;
            
            
            % add the current time when this function was called and the time when the file was created to logged_code_structure
            log_code_obj.logged_code_structure(curr_log_code_struct_ind).logged_time = datestr(now);   % can do other formats, such as datestr(now, 31), but easy to convert later using the datestr command
            curr_file_info = dir(full_file_name);
            log_code_obj.logged_code_structure(curr_log_code_struct_ind).last_modified_time = curr_file_info.date;
            
            
            % add the name of the file and the path to the file to the logged_code_structure
            [curr_path_string, curr_file_name, file_extension] = fileparts(full_file_name);
            log_code_obj.logged_code_structure(curr_log_code_struct_ind).code_filename = [curr_file_name file_extension];
            
            if isempty(curr_path_string)
                log_code_obj.logged_code_structure(curr_log_code_struct_ind).code_filepath = pwd;
            else
                log_code_obj.logged_code_structure(curr_log_code_struct_ind).code_filepath = curr_path_string;
            end
            
             
             % add the body of the code to the logged_code_structure
            f = fopen(full_file_name, 'r');
            if f < 0, error('unable to open "%s"', full_file_name); end
            curr_code_body = fread(f, '*char')';
            fclose(f);
            log_code_obj.logged_code_structure(curr_log_code_struct_ind).code_body = strrep(curr_code_body, sprintf('\r\n'), sprintf('\n'));    % replace \r\n new lines characters with \n new line characters  (not sure if this is really needed)

        end
        
          
    end   % end private methods
    
    
end

