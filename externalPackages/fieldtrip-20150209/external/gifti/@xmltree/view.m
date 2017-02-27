function view(tree)
% XMLTREE/VIEW View Method
% FORMAT view(tree)
% 
% tree   - XMLTree object
%__________________________________________________________________________
%
% Display an XML tree in a graphical interface
%__________________________________________________________________________
% Copyright (C) 2002-2011  http://www.artefact.tk/

% Guillaume Flandin
% $Id: view.m 8776 2013-11-14 09:04:48Z roboos $


%error(nargchk(1,1,nargin));

editor(tree);
