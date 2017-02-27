function [report_data] = report_objective(report_data, state)
%

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: report_objective.m 8776 2013-11-14 09:04:48Z roboos $

iteration = state.iteration;
component = state.component;
[params, s] = feval(state.denf.h, state.denf.params, state.s, state);
report_data.objective(component,iteration) = state.s * s';

% TODO: normalization

