function MCC = Restrict_Points(MCC,idx)
% INPUT: idx -- a list of indices (points)
% OUTPUT: MCC object 
%
% Create a list of illegal points (points not allowed in the cluster.
%
% cowen 2002 restrict results to points OTHER than these, regardless of convex hulls.
% modified ncst 02 May 03

MCC.ForbiddenPoints = unique([MCC.ForbiddenPoints(:); idx(:)]); 
MCC.recalc = 1;