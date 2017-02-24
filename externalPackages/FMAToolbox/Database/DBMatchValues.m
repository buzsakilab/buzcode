function [in1,in2,out1,out2,eid1,eid2] = DBMatchValues(eid1,eid2,rexp)

%DBMatchValues - Match values for two different sets of variables.
%
%  It is often necessary to match (associate) variables that were computed
%  independently, e.g. the phase precession slope and field size of hippocampal
%  place cells. Typically, such data are computed over repeated experiments
%  using batch processing tools (see <a href="StartBatch">StartBatch</a>) and stored in databases for
%  subsequent group analyses. Once retrieved using <a href="DBGetValues">DBGetValues</a>, they can be
%  matched using this function.
%
%  USAGE
%
%    [in1,in2,out1,out2,m1,m2] = DBMatchValues(eid1,eid2,rexp)
%
%    eid1,eid2      experiment IDs (see <a href="DBGetValues">DBGetValues</a>)
%    rexp           optional regular expression (see Example 2)
%
%  OUTPUT
%
%    in1            indices of elements in eid1 that are also in eid2
%    in2            indices of elements in eid2 that are also in eid1
%    out1           indices of elements in eid1 that are not in eid2
%    out2           indices of elements in eid2 that are not in eid1
%    m1             matches (extracted patterns) in eid1
%    m2             matches (extracted patterns) in eid2
%
%  EXAMPLE 1
%
%    In this example, place cells were recorded on successive days while the
%    animal explored a maze. Phase precession slopes and firing field sizes
%    were stored in a database using eids like '20120213-Maze-(1,2)' (where 1,2
%    corresponds to tetrode 1, cluster 2) and named 'PPSlope' and 'FieldSize',
%    respectively.
%
%    Get slopes and sizes:
%
%    [slopes,eid1] = DBGetValues('eid like "%Maze%" and name="PPSlope"');
%    [sizes,eid2] = DBGetValues('eid like "%Maze%" and name="FieldSize"');
%
%    Discard values that cannot be correlated (incomplete pairs), and reorder
%    so that each line in both variables corresponds to the same data:
%
%    [in1,in2] = DBMatchValues(eid1,eid2);
%    slopes = slopes(in1);
%    sizes = sizes(in2);
%
%  EXAMPLE 2
%
%    In this example, cells were recorded during both wake an sleep, and their
%    average firing rates should be compared across behavioral conditions.
%    Data were stored in a database using eids like '20120213-Maze-(1,2)' or
%    '20130214-Sleep-(1,2)' and named 'MeanRate'. Matching them is trickier
%    because their eids are not pairwise identical. Here we need to extract
%    the relevant portions of the eids, i.e. discard 'Sleep' or 'Maze'.
%
%    Get mean rates:
%
%    [maze,eid1] = DBGetValues('eid like "%Maze%" and name="MeanRate"');
%    [sleep,eid2] = DBGetValues('eid like "%Sleep%" and name="MeanRate"');
%
%    Discard values that cannot be correlated (incomplete pairs), and reorder
%    so that each line in both variables corresponds to the same data:
%
%    [in1,in2] = DBMatchValues(eid1,eid2,'([0-9]{8}).*([0-9]*,[0-9]*)');
%    maze = maze(in1);
%    sleep = sleep(in2);
%
%    The regular expression includes two tokens (indicated by the parentheses),
%    namely [0-9]{8} (eight successive occurrences of a digit between 0 and 9)
%    and [0-9]*,[0-9]* (any number of digits, a comma, any number of digits).
%
%  SEE
%
%    See also DBGetValues, DBGetVariables, DBAddVariable.
%


% Copyright (C) 2013 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Transform eids if necessary
if nargin >= 3,
	eid1 = DoRegexp(eid1,rexp);
	eid2 = DoRegexp(eid2,rexp);
else
	disp([ 'Example eid: ' eid1{1}]);
	disp([ 'Example eid: ' eid2{1}]);
end
% Make sure at least one of the lists contains unique values (keys)
key1 = length(eid1) == length(unique(eid1));
key2 = length(eid2) == length(unique(eid2));
if ~key1 && ~key2,
	error('Neither list contains unique values.');
end

% Check which eids in list 1 are in list 2, making sure empty eids always fail the test
[~,in2] = ismember(eid1,eid2);
empty = cellfun('isempty',eid1);
in2(empty) = 0;
% The code above returns locations (items found) intermixed with zeros (items not found)
% (for details, see help for 'ismember'). Split this information into separate variables.
out1 = find(in2==0);
in2 = in2(in2~=0);
% If values in list 2 are not unique, duplicate entries in list 1 where necessary
if ~key2,
	in2 = find(ismember(eid2,eid2(in2)));
end

% Check which eids in list 2 are in list 1, making sure empty eids always fail the test
[~,in1] = ismember(eid2,eid1);
empty = cellfun('isempty',eid2);
in1(empty) = 0;
% The code above returns locations (items found) intermixed with zeros (items not found)
% (for details, see help for 'ismember'). Split this information into separate variables.
out2 = find(in1==0);
in1 = in1(in1~=0);
% If values in list 1 are not unique, duplicate entries in list 2 where necessary
if ~key1,
	in1 = find(ismember(eid1,eid1(in1)));
end

% Helper function: transform all eids, i.e. extract tokens using the regular expression then concatenate
% the tokens (eids that do not match are set to the empty string '')

function r = DoRegexp(s,rexp)

r = cellfun(@(x) regexp(x,rexp,'tokens'),s,'uniformoutput',false);
empty = cellfun('isempty',r);
r(empty) = {''};
r(~empty) = cellfun(@(x) horzcat(x{1}{:}),r(~empty),'uniformoutput',false);

disp([ 'Example pattern extraction: ' s{1} ' -> ' r{1}]);