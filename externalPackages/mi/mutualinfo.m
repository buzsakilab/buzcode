function h = mutualinfo(vec1,vec2)
%=========================================================
%
%This is a prog in the MutualInfo 0.9 package written by 
% Hanchuan Peng.
%
%Disclaimer: The author of program is Hanchuan Peng
%      at <penghanchuan@yahoo.com> and <phc@cbmv.jhu.edu>.
%
%https://www.mathworks.com/matlabcentral/fileexchange/14888-mutual-information-computation?s_tid=FX_rc1_behav
%The CopyRight is reserved by the author.
%
%Last modification: April/19/2002
%
%========================================================
%
% h = mutualinfo(vec1,vec2)
% calculate the mutual information of two vectors
% By Hanchuan Peng, April/2002
%

%% Fix for inf/nan bug
infnan = isinf(vec1) | isinf(vec2) | isnan(vec1) | isnan(vec2);
vec1(infnan) = [];
vec2(infnan) = [];

%% FIx single bug
if isa(vec1,'single')
    vec1 = double(vec1);
end
if isa(vec2,'single')
    vec2 = double(vec2);
end

if isempty(vec1) || isempty(vec2)
   h = nan;
   return
end
%%
[p12, p1, p2] = estpab(vec1,vec2);
h = estmutualinfo(p12,p1,p2);

