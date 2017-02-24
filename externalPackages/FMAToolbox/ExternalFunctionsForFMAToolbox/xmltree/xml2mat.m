function varargout = xml2mat(xmlfile, matfile)
%XML2MAT Convert an XML-file into a MAT-file
%  XML2MAT(XMLFILE, MATFILE) converts an XML file or an XML string
%  (as saved by MAT2XML) into a MAT-file MATFILE.
%  S = XML2MAT(XMLFILE) returns the content of XMLFILE in variable S.
%  S is a struct containing fields matching the variables retrieved.
%
%  See also LOADXML, MAT2XML, XMLTREE, LOAD.

%  Copyright 2003 Guillaume Flandin. 
%  Revision: 1.0 $  $Date: 2003/07/01 21:43 $

% Bugs to handle:
%  o xml_parser replaces tabs and successive spaces into one even
%    in CData tags (=> pbs for char variables)
%  o still a bug with entities ?
%  o integer variables not handled
%  o proper handle of home made classes
%  o deal with precision loss with double variables (format)
%  o special handling of cellstr and struct array ?
%  o does it work with all Matlab releases ?
%  o if xmlfile is "empty" and nargout = 1 then error
%  o Use of 'eval' may reinitialize variables used in this function

error(nargchk(1,2,nargin));

s = loadxml(xmlfile);

if nargout == 1 | nargin == 1
	varargout{1} = s;
end

if nargin == 2
	names = fieldnames(s);
	flagfirstvar = 1;
	for i=1:length(names)
		% TODO % Very Dangerous !!!
		eval([names{i} ' = s.' names{i} ';']);
		if flagfirstvar
			save(matfile,names{i});
			flagfirstvar = 0;
		else
			save(matfile,names{i},'-APPEND');
		end
	end
end
