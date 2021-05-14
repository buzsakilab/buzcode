function [ fields ] = bz_getPlaceFields(varargin)
% USAGE
%
%
% INPUTS
% 
%   ratemap   MxNxD matrix where M is the number of cells, N is the number
%             of trials, and D is the number of spatial bins
%
%
% OUTPUTS
%
%   fields struct with field data for each cell
%
%
%
% HELP
% This function tries to identify place fields based on a set of heuristics
% 
% written by david tingley, 2017

p = inputParser;
addRequired(p,'ratemap',1250,@isnumeric)
addParameter(p,'minPeakRate',5,@isnumeric)
addParameter(p,'minFieldWidth',8,@isnumeric)
addParameter(p,'maxFieldWidth',100,@isnumeric)

% addParameter(p,'minPeakRate',3,@isnumeric)

ratemap = p.Results.ratemap;
minPeakRate = p.Results.minPeakRate;
minFieldWidth = p.Results.minFieldWidth;
maxFieldWidth = p.Results.maxFieldWidth;

%% find peaks in avg firing rates above 3 hz

meanRates = squeeze(mean(ratemap,2));


for i=1:size(meanRates,1)
    [pks locs w] = findpeaks(fastrms(meanRates(i,:),5),'minpeakheight',minPeakRate,'MinPeakWidth',minFieldWidth);
    exclude=[];
    for j=1:length(locs)-1
       if min(meanRates(i,locs(j):locs(j+1))) > ((pks(j)+pks(j+1))./2) * .1
           % exclude fields without a 70% decrease in rate between peaks
           if pks(j) > pks(j+1)
               exclude = [exclude;j+1];
           elseif pks(j) < pks(j+1)
               exclude = [exclude;j];
           end
       end
    end   
    pks(exclude) = [];
    locs(exclude)=[];
   
    for j=1:length(locs)
        
        Map_Field = meanRates(i,:) > pks(j) * .1;
        
        start = locs(j);
        stop = locs(j);
        while Map_Field(start) == 1
            start = start-1;
        end
        while Map_Field(stop) == 1
            stop = stop+1;
        end
        if stop - start > minFieldWidth && stop - start < maxFieldWidth
            fields{i}{j}.start = start;
            fields{i}{j}.stop = stop;
            fields{i}{j}.width = stop - start;
            fields{i}{j}.peakFR = pks(j);
            fields{i}{j}.peakLoc = locs(j);
        end
    end
end






end

