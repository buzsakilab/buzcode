function [ paddeddata ] = NanPadJumps( timestamps,data,minjump )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
jumpstarts = find(diff(timestamps)>minjump);
jumpstops = jumpstarts+1;

paddeddata = data;
paddeddata([jumpstarts jumpstops]) = nan;

%%
% figure
% plot(timestamps,paddeddata)
% hold on
% plot(timestamps(isnan(paddeddata)),...
%     ones(size(timestamps(isnan(paddeddata)))),'+')
end

