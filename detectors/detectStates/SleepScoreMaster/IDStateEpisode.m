function [ STATEints ] = IDStateEpisode(STATEints,maxintdur,minSTATEdur)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%INPUTS
%   STATEints
%   maxintdur   maximum interruption
%   minSTATEdur minimum state duration


%%
%SWSints = stateintervals{2};

[STATEints] = MergeSeparatedInts_ss(STATEints,maxintdur);
STATElengths = STATEints(:,2)-STATEints(:,1);
STATEints = STATEints(STATElengths>=minSTATEdur,:);




        
end

