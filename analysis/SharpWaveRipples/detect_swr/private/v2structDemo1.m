function [] = v2structDemo1()
%  V2STRUCTDEMO1 is a demo function showing some examples for using the function V2STRUCT 
%
%% Last update:
%    30.06.2011, Adi N.

%%
fprintf(['\nRunning demo function with examples for using the function vs2truct:\n',...
        'Press any key between breaks.\n']);
pause()
fprintf(['\nInitialize demo variables:\n',...
         'x = zeros(3); x2 = ones(3); y = ''Testing123''; z = cell(2,3);\n',...
         'fieldNames1 = {''fieldNames'',''x'',''y'',''z''};\n',...
         'fieldNames2 = {''fieldNames'',''a'',''b'',''c''};\n',...
         'fieldNames3 = {''fieldNames'',''x''};\n',...
         'nameOfStruct2Update = ''S'';\n']);
pause()
x = zeros(3); x2 = ones(3); y = 'Testing123'; z = cell(2,3);
fieldNames1 = {'fieldNames','x','y','z'};
fieldNames2 = {'fieldNames','a','b','c'};
fieldNames3 = {'fieldNames','x'};
nameOfStruct2Update = 'S';

fprintf(['\nExamples for packing and unpacking:\n',...
         '\nPacking:\n']);
pause()
fprintf('\nExample 1: Variables x,y and z are packed to structure S\n');
fprintf('S = v2struct(x,y,z) ');
pause()
S = v2struct(x,y,z) 

pause()
fprintf(['Example 2: Field names defined according to the cell array fieldNames.\n',...
         'variables with the names in fieldNames must exist.\n']);
clear S
fprintf('S = v2struct(fieldNames1)');
pause()
S = v2struct(fieldNames1)

pause()
fprintf('Example 3: Same as the first but arguments are passed explicitly\n');
clear S
fprintf('S = v2struct(zeros(3),''Testing123'',cell(2,3),fieldNames1)');
pause()
S = v2struct(zeros(3),'Testing123',cell(2,3),fieldNames1)

pause()
fprintf(['Example 4: Field names defined by content of fieldNames2 while the values are\n',...
         'set according to the passed arguments. In this case the structure S returned would\n'...
         'be: S.a=x, S.b=y, S.c=z\n']);
clear S
fprintf('S = v2struct(x,y,z,fieldNames2)');
pause()
S = v2struct(x,y,z,fieldNames2)

pause()
fprintf(['Example 5: Update structure S. The fields that would be updated are according to content\n',...
         'of fieldNames3. Note that you must pass a variable with the name\n',...
         '''nameOfStruct2Update'' placed before ''fieldNames3''. This variable should\n',...
         'contain the name of the structure you want to update as a string. Also note\n',...
         'that if you set an output structure name which is different than the one stated\n'...
         'in nameOfStruct2Update a new structure would be created and the structure that\n'...
         'was meant to be updated would not get updated.']);
fprintf(['S.oldField = ''field to be saved for future use''\n',...
        'S = v2struct(x2,nameOfStruct2Update,fieldNames3)\n']);
pause()
clear S
S.oldField = 'field to be saved for future use'
S = v2struct(x2,nameOfStruct2Update,fieldNames3)

pause()
fprintf('Example 6: Pack all variables in caller workspace.\n');
fprintf('S = v2struct');
pause()
S = v2struct

%%
pause()
fprintf(['\nExamples for Unpacking:\n']);
pause()

fprintf(['\nFisrt clear variables from ''packing'' and initialize variables for\n',...
         '''unpacking'':\n',...
         'clear all\n',...
         'S.x = zeros(3); S.y = ''Testing123''; S.z = cell(2,3);\n',...
         'fieldNames3 = {''fieldNames'',''x'',''z''};\n'])

clear all
S.x = zeros(3); S.y = 'Testing123'; S.z = cell(2,3); fieldNames3 = {'fieldNames','x','z'};

pause()
fprintf(['\nExample 1: Create or overwrite variables x, y, z in the caller environment\n',...
         'with the contents of the corresponding named fields.\n']);
pause()
fprintf(['v2struct(S)\n']);
pause()

v2struct(S)

fprintf(['\nExample 2: Assign the contents of the fields of the scalar structure S to\n',...
         'the variables a,b,c rather than overwriting variables in the caller. If there\n',...
         'are fewer output variables than there are fields in S, the remaining fields\n',...
         'are not extracted.\n']);
pause()
fprintf(['[a,b,c] = v2struct(S)\n']);

pause()

[a,b,c] = v2struct(S)

fprintf(['\nExample 3: Create or overwrite variables x and z in the caller with the\n',...
         'contents of the corresponding named fields.\n']);
pause()
fprintf(['v2struct(S,fieldNames3)\n']);
pause()
v2struct(S,fieldNames3)

fprintf(['\nExample 4: Assign the contents of the fields ''x'' and ''z'' defined by\n',...
         'fieldNames3 of the scalar structure S to the variables a and b rather than\n',...
         'overwriting variables in the caller. If there are fewer output variables than\n',...
         'there are fields in S, the remaining fields are not extracted.\n']);
pause()
fprintf(['[a,b] = v2struct(S,fieldNames3)\n']);
pause()
[a,b] = v2struct(S,fieldNames3)

fprintf(['\nExample 5: Unpack variable ''z'' only without overwriting variable ''x''. Note\n',...
         'the addition of the field named ''avoidOverWrite'' to the structure to be\n',...
         'unpacked. This is mandatory in order to make this functionality work. The\n',...
         'contents of this field can be anything, it does not matter.\n']);

pause()
fprintf(['S.avoidOverWrite = ''foo(contents does not matter)'';\n',...
         'x = ''do not overwrite me'';\n',...
         'clear y\n']);
pause()
S.avoidOverWrite = 'foo(contents does not matter)';
x = 'do not overwrite me';
v2struct(S)
S

fprintf('\nend of v2structDemo1\n');
