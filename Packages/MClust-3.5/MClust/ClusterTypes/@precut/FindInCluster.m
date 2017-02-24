function [f, MCC] = FindInCluster(MCC)

% [f, MCC] = FindInCluster(MCC, varargin)
% Find in cluster from already cut data
% added sort ADR 16 may 2002
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.

f = MCC.myPoints;
f = sort(f);