function [report_data] = report_w(report_data, state)
% Reporting function to record w for each iteration
%

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: report_w.m 8776 2013-11-14 09:04:48Z roboos $

report_data.w(:, state.iteration) = state.w;
