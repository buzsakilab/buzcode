%% v2struct
% v2struct Pack/Unpack Variables to/from a scalar structure.
function varargout = v2struct(varargin)

%% Description
%    v2struct has dual functionality in packing & unpacking variables into structures and
%    vice versa, according to the syntax and inputs.
%
%    Function features:
%      * Pack variables to structure with enhanced field naming
%      * Pack and update variables in existing structure
%      * Unpack variables from structure with enhanced variable naming
%      * Unpack only specific fields in a structure to variables
%      * Unpack without over writing existing variables in work space
%
%    In addition to the obvious usage, this function could by highly useful for example in
%    working with a function with multiple inputs. Packing variables before the call to
%    the function, and unpacking it in the beginning of the function will make the function
%    call shorter, more readable, and you would not have to worry about arguments order any
%    more. Moreover you could leave the function as it is and you could pass same inputs to
%    multiple functions, each of which will use its designated arguments placed in the
%    structure.
%
%% Syntax
%    Pack
%      S = v2struct
%      S = v2struct(x,y,z,...)
%      S = v2struct(fieldNames)
%      S = v2struct(A,B,C,..., fieldNames)
%      S = v2struct(x,..., nameOfStruct2Update, fieldNames)
%      v2struct
%      v2struct(x,y,z,...)
%      v2struct(fieldNames)
%      v2struct(A,B,C,..., fieldNames)
%      v2struct(x,..., nameOfStruct2Update, fieldNames)
%
%    Unpack
%      v2struct(S)
%      [a,b,c,...] = v2struct(S)
%      v2struct(S,fieldNames)
%      [a,b,c,...] = v2struct(S,fieldNames)
%
%% Inputs & Outputs
%    Pack - inputs
%      x,y,z,...           - any variable to pack. can be replaced by fieldNames below.
%      nameOfStruct2Update - optional, name of structure to update if desired.
%      fieldNames          - optional, cell array of strings, which must include a cell
%                            with the string 'fieldNames' and must be the last input.
%    Pack - outputs 
%      S - the packed structure. If there is no output argument then a structure named
%          Sv2struct would be created in the caller workspace.
%
%    Unpack - inputs
%      S          - name of structure to be unpacked.
%      fieldNames - optional, cell array of strings, which must include a cell with the
%                   string 'fieldNames' and must be the last input.
%    Unpack - outputs          
%      a,b,c,... - variables upacked from the structure.
%                  if there are no output arguments then variables would be created in
%                  the caller workspace with naming according to name of inputs.
%
%% Examples
%  % see 'Usage example' section below for convenient presentation of these examples.
%
%    % NOTE: whenever using filedNames cell array please note the following
%    %  1. fieldNames cell array must include a cell with the string 'fieldNames'
%    %  2. fieldNames cell array input must be the last input.
%
%  % Pack
%      x = zeros(3); x2 = ones(3); y = 'Testing123'; z = cell(2,3);
%      fieldNames1 = {'fieldNames', 'x', 'y', 'z'};
%      fieldNames2 = {'fieldNames', 'a', 'b', 'c'};
%      fieldNames3 = {'fieldNames', 'x'};
%      nameOfStruct2Update = 'S';
%
%       % The four examples below return structure S with same values however the
%       % structure's field names are defined differently in every syntax. 
%      % Example 1.
%      % structure field names defined by variables names.
%       S = v2struct(x,y,z) 
%      % Example 2.
%      % structure field names defined according to the cell array fieldNames. 
%       % NOTE: variables with the names in fieldNames1 must exist in the caller workspace.
%       S = v2struct(fieldNames1) 
%      % Example 3.
%      % same as #1. but arguments are passed explicitly
%       S = v2struct(zeros(3), 'Testing123', cell(2,3), fieldNames1) 
%      % Example 4.
%      % field names defined by content of fieldNames2 while
%      % the values are set according to the passed arguments. In this case the structure
%      % S returned would be: S.a=x, S.b=y, S.c=z
%       S = v2struct(x,y,z, fieldNames2) 
%
%      % Example 5.
%      % update structure S. The fields that would be updated are according to content
%      % of fieldNames3. Note that you must pass a variable with the name
%      % 'nameOfStruct2Update' placed before 'fieldNames3'. This variable should contain
%      % the name of the structure you want to update as a string. Also note that if you
%      % set an output structure name which is different than the one stated in
%      % nameOfStruct2Update a new structure would be created and the structure that was
%      % meant to be updated would not get updated.
%       S.oldField = 'field to be saved for future use'
%       S = v2struct(x2, nameOfStruct2Update, fieldNames3)
%
%      % Example 6.
%      % pack all variables in caller workspace. Call without input arguments.
%        S = v2struct
%
%      % The following examples return the same results as the examples above but the
%      % structure would be returned with the default name 'Sv2struct'. Be cautious as
%      % this might lead to overriding of arguments.
%      % Example 7.
%       v2struct(x,y,z)
%      % Example 8.
%       v2struct(fieldNames1)
%      % Example 9.
%       v2struct(zeros(3), 'Testing123', cell(2,3), fieldNames1)
%      % Example 10.
%       v2struct(x,y,z, fieldNames2)
%      % Example 11.
%       S.oldField = 'field to be saved for future use'
%       v2struct(x2, nameOfStruct2Update, fieldNames3)
%      % Example 12.
%       v2struct
%
%  % Unpack
%      clear S x x2 y z fieldNames1 fieldNames2 fieldNames3 nameOfStruct2Update
%      S.x = zeros(3); S.y = 'Testing123'; S.z = cell(2,3);
%      fieldNames3 = {'fieldNames','x','z'};
%
%      % Example 1.
%      % This example creates or overwrites variables x, y, z in the caller with the
%      % contents of the corresponding named fields.
%       v2struct(S)
%
%      % Example 2.
%      % This example assigns the contents of the fields of the scalar structure
%      % S to the variables a,b,c rather than overwriting variables in the caller. If
%      % there are fewer output variables than there are fields in S, the remaining fields
%      % are not extracted.
%       [a,b,c] = v2struct(S)
%
%      % Example 3.
%      % This example creates or overwrites variables x and z in the caller with the
%      % contents of the corresponding named fields.
%       v2struct(S, fieldNames3)
%
%      % Example 4.
%      % This example assigns the contents of the fields 'x' and 'z' defined by
%      % fieldNames3 of the scalar structure S to the variables a and b rather than
%      % overwriting variables in the caller. If there are fewer output variables than
%      % there are fields in S, the remaining fields are not extracted.
%       [a,b] = v2struct(S, fieldNames3)
%
%       % This example unpacks variables 'y' and 'z' only without overwriting variable 'x'. 
%       % NOTE the addition of the field named 'avoidOverWrite' to the structure to be
%       % unpacked. This is mandatory in order to make this functionality work. The
%       % contents of this field can be anything, it does not matter. 
%      S.avoidOverWrite = '';
%      x = 'do not overwrite me';
%      v2struct(S)
%
%% Usage example (includes sub-functions)
%    1. run attached v2structDemo1.m file for on screen presentation of examples.
%    2. run attached v2structDemo2.m file and read comments in file for a suggestion of
%       how to use v2struct in managing input to other functions with improved usability.
%
%% Revision history
%    2011-05-19, Adi N., Creation
%    2011-05-29, Adi N., Update structure added, some documentation and demo function changes
%    2011-06-02, Adi N., Fixed updating structure functionality
%    2011-06-05, Adi N., Added functionality: avoid overwritring existing variables, added
%                        unpacking examples to demo1 .m file.
%    2011-06-30, Adi N., fieldNames usage corrected, now must include a specific string to
%                        be triggered. Documentation enhanced. Code tweaked.
%    2011-07-14, Adi N., Fixed bug in packing with variables only
%    2011-08-14, Adi N., Clarified warning and error when packing/unpacking with
%                        fieldNames.
%    2011-09-12, Adi N., Added easy packing of all variables in caller workspace (thanks 
%                        to Vesa Lehtinen for the suggestion), fixed bug in warning
%                        handling in packing case, edited comments.
%
%    Inspired by the function: mmv2truct - D.C. Hanselman, University of Maine, Orono, ME
%    04469 4/28/99, 9/29/99, renamed 10/19/99 Mastering MATLAB 5, Prentice Hall,
%    ISBN 0-13-858366-8


% parse input for field names
if isempty(varargin)
   gotCellArrayOfStrings = false;
   toUnpackRegular = false;
   toUnpackFieldNames = false;
   gotFieldNames = false;
else
   gotCellArrayOfStrings = iscellstr(varargin{end});
   
   toUnpackRegular = (nargin == 1) && isstruct(varargin{1});
   if toUnpackRegular
      fieldNames = fieldnames(varargin{1})';
      nFields = length(fieldNames);
   end
   
   gotFieldNames = gotCellArrayOfStrings & any(strcmpi(varargin{end},'fieldNames'));
   if gotFieldNames
      fieldNamesRaw = varargin{end};
      % indices of cells with actual field names, leaving out the index to 'fieldNames' cell.
      indFieldNames = ~strcmpi(fieldNamesRaw,'fieldNames');
      fieldNames = fieldNamesRaw(indFieldNames);
      nFields = length(fieldNames);
   end
   toUnpackFieldNames = (nargin == 2) && isstruct(varargin{1}) && gotFieldNames;
end


% Unpack
if toUnpackRegular || toUnpackFieldNames 
   
   struct = varargin{1};
   assert(isequal(length(struct),1) , 'Single input nust be a scalar structure.');
   CallerWS = evalin('caller','whos'); % arguments in caller work space
   
   % update fieldNames according to 'avoidOverWrite' flag field.
   if isfield(struct,'avoidOverWrite')
      indFieldNames = ~ismember(fieldNames,{CallerWS(:).name,'avoidOverWrite'});
      fieldNames = fieldNames(indFieldNames);
      nFields = length(fieldNames);
   end
   
   if toUnpackRegular % Unpack with regular fields order
      if nargout == 0 % assign in caller
         for iField = 1:nFields
            assignin('caller',fieldNames{iField},struct.(fieldNames{iField}));
         end
      else % dump into variables
         for iField = 1:nargout
            varargout{iField} = struct.(fieldNames{iField});
         end
      end
      
   elseif toUnpackFieldNames % Unpack with fields according to fieldNames
      if nargout == 0 % assign in caller, by comparing fields to fieldNames
         for iField = 1:nFields
            assignin('caller',fieldNames{iField},struct.(fieldNames{iField}));
         end
      else % dump into variables
         assert( isequal(nFields, nargout) , ['Number of output arguments',...
            ' does not match number of field names in cell array']);
         for iField = 1:nFields
            varargout{iField} = struct.(fieldNames{iField});
         end
      end
   end
   
% Pack   
else
   % build cell array of input names
   CallerWS = evalin('caller','whos');
   inputNames = cell(1,nargin);
   for iArgin = 1:nargin
      inputNames{iArgin} = inputname(iArgin);
   end
   nInputs = length(inputNames);
   
      % look for 'nameOfStruct2Update' variable and get the structure name
   if ~any(strcmpi(inputNames,'nameOfStruct2Update')) % no nameOfStruct2Update
      nameStructArgFound = false;
      validVarargin = varargin;
   else % nameOfStruct2Update found
      nameStructArgFound = true;
      nameStructArgLoc = strcmp(inputNames,'nameOfStruct2Update');
      nameOfStruct2Update = varargin{nameStructArgLoc};
      % valid varargin with just the inputs to pack and fieldNames if exists
      validVarargin = varargin(~strcmpi(inputNames,'nameOfStruct2Update'));
      % valid inputNames with just the inputs name to pack and fieldNames if exists
      inputNames = inputNames(~strcmpi(inputNames,'nameOfStruct2Update'));
      nInputs = length(inputNames);
      % copy structure from caller workspace to enable its updating
      if ismember(nameOfStruct2Update,{CallerWS(:).name}) % verify existance
         S = evalin('caller',nameOfStruct2Update);
      else
         error(['Bad input. Structure named ''',nameOfStruct2Update,...
            ''' was not found in workspace'])
      end
   end
   
   % when there is no input or the input is only variables and perhaps
   % also nameOfStruct2Update   
   if ~gotFieldNames
      % no input, pack all of variables in caller workspace
      if isequal(nInputs, 0)
         for iVar = 1:length(CallerWS)
            S.(CallerWS(iVar).name) = evalin('caller',CallerWS(iVar).name);
         end
         % got input, check input names and pack
      else
         for iInput = 1:nInputs
            if gotCellArrayOfStrings % called with a cell array of strings
               errMsg = sprintf(['Bad input in cell array of strings.'...
                           '\nIf you want to pack (or unpack) using this cell array as'...
                           ' designated names'...
                           '\nof the structure''s fields, add a cell with the string'...
                           ' ''fieldNames'' to it.']);
            else
               errMsg = sprintf(['Bad input in argument no. ', int2str(iArgin),...
                                 ' - explicit argument.\n'...
                          'Explicit arguments can only be called along with a matching'...
                          '\n''fieldNames'' cell array of strings.']);
            end
            assert( ~isempty(inputNames{iInput}), errMsg);
            S.(inputNames{iInput}) = validVarargin{iInput};
         end
         
         % issue warning for possible wrong usage when packing with an input of cell array of
         % strings as the last input without it containing the string 'fieldNames'.
         if gotCellArrayOfStrings
            name = inputNames{end};
            % input contains structure and a cell array of strings
            if (nargin == 2) && isstruct(varargin{1})
               msgStr = [inputNames{1},''' and ''',inputNames{2},''' were'];
               % input contains any arguments with an implicit cell array of strings
            else
               msgStr = [name, ''' was'];
            end
            warnMsg = ['V2STRUCT - ''%s packed in the structure.'...
               '\nTo avoid this warning do not put ''%s'' as last v2struct input.'...
               '\nIf you want to pack (or unpack) using ''%s'' as designated names'...
               ' of the'...
               '\nstructure''s fields, add a cell with the string ''fieldNames'' to'...
               ' ''%s''.'];
            fprintf('\n')
            warning('MATLAB:V2STRUCT:cellArrayOfStringNotFieldNames',warnMsg,msgStr,...
                     name,name,name)
         end
      end
   % fieldNames cell array exists in input
   elseif gotFieldNames
      nVarToPack = length(varargin)-1-double(nameStructArgFound);
      if nVarToPack == 0 % no variables to pack
         for iField = 1:nFields
            S.(fieldNames{iField}) = evalin('caller',fieldNames{iField});
         end
         
         % else - variables to pack exist
         % check for correct number of fields vs. variables to pack
      elseif ~isequal(nFields,nVarToPack)
         error(['Bad input. Number of strings in fieldNames does not match',...
                'number of input arguments for packing.'])
      else
         for iField = 1:nFields
            S.(fieldNames{iField}) = validVarargin{iField};
         end
      end
      
   end % if ~gotFieldNames

if nargout == 0
   assignin( 'caller', 'Sv2struct',S );
else
   varargout{1} = S;
end

end % if nargin
