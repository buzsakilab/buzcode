% Several examples for using StructDlg

%%
% Note: The values the user enters to the dialog box for numeric fields are being 'eval'-uated,
% so inputs such as (sin(0.34)+2)^my_func(34.2) are possible.
% The field's numeric limits can be displayed by placing the mouse over the field.
%
%
% Free Strings: Limited to one line string (no string arrays). If no default is required, use the  
%                     empty string (''), otherwise the field will be treated as numeric.

clear S
S.center_frequency = { 2000 'Hz' [30 50000] }; %  default of 2000, allowed range:[30 50000].
S.my_parameter     = { 43 }; %  default of 43, no limits.
S.size_of_matrix   = { [3 2] '' [1 Inf] }; %  default of [4 12], allowed range:[1 Inf].
S.Weight_Matrix = [0.3 0.4; 3.4 9.1; 10 0.4];
S.User_name = {'' ''};
P = StructDlg(S,'My title')

%%
% Radio buttons: The definition includes a string that contains the options
% separated by '|'. The option that the user chooses will become the value 
%                     of that field in P. A Default option can be specified by enclosing it with curly brackets. 
%                     The chosen option is converted to a numeric value if possible.
clear S;
S.colormap = {'hsv|{gray}|hot|bone|pink'};
S.Sampling_rate = {'11500|22000|{44000}', 'Hz'};
P = StructDlg(S,'Radio buttons')

%%
% Pop-up menu: The definition includes a Cell-array of single strings. 
%                     The chosen option will be the value of that field in P. A default option can be 
%                     specified by enclosing it with curly brackets. The chosen string is converted
%                     to a numeric value if possible.
%
% Check box: The definition must be in the form: {'0' '1'}. 
%                A default value may be specified by the curly brackets.

clear S;
S.Sampling_rate = { {'12207' '24414' '48828' '{97656}'} , 'Hz'};
S.Filter_type = { {'bartlett' '{chebwin}' 'hamming'} };
S.use_filter = { {'0','{1}'} };
P = StructDlg(S,'Popup menus and check box')
         
%%
% File Name Dialogs: The definition includes a cell of one string that must start with 'uigetfile' or
%                     'uiputfile'. The getfile commands may be followed by a default file name or
%                     a search filter enclosed by brackets. The user will be able to specify the file
%                     name directly by typing it, or to push a small pushbutton that will pop-up 
%                     Matlab's uigetfile or uipufile.
clear S
S.parameters_file = { {'uigetfile(''c:\temp\*.mat'')'} };
P = StructDlg(S,'File name dialogs');
%%
% Sub-Structure: S may contain substrucutres of the same format. The user will be able to push a
%                     push-button that will call 'StructDlg' recursively for the sub-structure.
%                     The current values of the sub-structure can be viewed by placing the mouse over
%                     the push-button.

%%
% Field Self-Reference: (For numeric fields only) A numeric field may include a reference to the value of 
%                     another numeric fields. This is done using a reserved word 'this' to refer to the
%                     structure.
clear S
S.Window_size       = { 512 'samples' [10 1000] };
S.Overlap           = {'this.Window_size / 2' '' [0 Inf]}; %  Note that a non-empty limits indicator is needed in order to indicate that this is a numeric field.
P = StructDlg(S,'Self reference fields');
%%
% Notes: The value of the dependent field will be automatically updated when the values
%                     of the referenced fields changes. The automatically changed value will blink twice
%                     to alert the user. The user is then able to undo the automatic change using the 
%                     mouse's right-click.
%                     It is not possile for a sub-structure field to reference fields in other sub-structures,
%                     or in other structure levels. It is possible for a field to reference fields in 
%                     sub-structures of lower levels; however, this is highly not recommended. 

%     'title' is the title of the dialog-box. 'dflt' is a structure which contain default values for P. 
%     These values will override  the default values specified in struct_def. 
%     Note that the dflt should contain the default values ONLY, not the entire 
%     definition as specified in S. This is useful for enabling user-specific defaults that override the 
%     'factory defaults'.

%%
%    Final Example: 
clear S;
S.Title       = { ['Anaysis results as of ' date] ''}; 
S.Input_file  = { {'uigetfile(''c:\*.xls'')'} };
S.Output_file = { {'uiputfile(''c:\*.txt'')'} };
S.Replace_NaNs_with_zeros = { {'{0}' '1'} };
S.Analysis_parameters.Method = { {'Linear' '{exponential}' 'logistic'} };
S.Analysis_parameters.Window_size = { 50 'samples' [5 100] };
S.Analysis_parameters.Theta = { 'min(this.Window_size, 10)' '' [0 100] };
P = StructDlg(S,'Sub Menus');
% after a while, we can ask again for user input.
% This time, the default values will be the ones the user specified previously.
% new_P = StructDlg(S,'Populattion analysis',P);
%

% clear S
% 
% S.Input_file  = { {'uigetfile(''c:\Images\*.jpg'')'} };
% S.Filter      = { {'{None}' 'Blur' 'Sharpen' 'Edge'} };
% S.BWconvert   = { {'{0}','1'} 'Convert to B/W'};  % default is NO
% S.quality     = { 80 'Output Quality' [1 100]};   % default of 80, allowed range: 1-100
% S.comment     = {''};
% S.Output_file = { {'uiputfile(''c:\Processed_Images\*.jpg'')'} };
% 
% P = StructDlg(S,'Parameters');
% disp(P)

clear S;
S.Name        = { 'Untitled' };
S.Width       = { 256 '* (pixels)' [1 1024] };             % default of 256, allowed range: 1-1024
S.Height      = { 256 '* (pixels)' [1 1024] };
S.Mode        = { {'Bitmap' 'Grayscale' '{RGB color}'} };  % default is 'RGB color'
S.Contents    = { '{White}|Background Color|Transperent'}; % default is 'White'
S.Use_alpha   = { {'{0}','1'} 'Use Alpha Channel' };       % default is NO
S.Output_file = { {'uiputfile(''c:\Images\*.jpg'')'} };

P = StructDlg(S,'New Image')
disp(P)