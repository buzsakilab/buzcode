function [ fields ] = bz_getPlaceFields1D(varargin)
% USAGE
%
%
% INPUTS
% 
%   ratemap             MxNxD matrix where M is the number of cells, N is the number
%                       of trials, and D is the number of spatial bins
%   minPeakRate         minimum rate for peak of field [default: 1]
%   minFieldWidth       minimum width of field [default: 5]
%   maxFieldWidth       maximum width of field [default: 100]
%   percentThreshold    percent change between peak rate and start/stop of field
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
addRequired(p,'ratemap',@isnumeric)
addParameter(p,'minPeakRate',1,@isnumeric)
addParameter(p,'minFieldWidth',5,@isnumeric)
addParameter(p,'maxFieldWidth',100,@isnumeric)
addParameter(p,'percentThreshold',.01,@isnumeric)

% addParameter(p,'minPeakRate',3,@isnumeric)
parse(p,varargin{:})

ratemap = p.Results.ratemap;
minPeakRate = p.Results.minPeakRate;
minFieldWidth = p.Results.minFieldWidth;
maxFieldWidth = p.Results.maxFieldWidth;
threshold = p.Results.percentThreshold; % change between peak rate and start/stop of field

%% find peaks in avg firing rates above 3 hz

meanRates = squeeze(mean(ratemap,2));

stdRates = squeeze(std(ratemap,[],2));
warning off  % findpeaks.m throws warnings if peak isn't found...

for i=1:size(meanRates,1)
    fields{i} = [];
    [pks locs w] = findpeaks(fastrms(meanRates(i,:),5),'minpeakheight',minPeakRate,'MinPeakWidth',minFieldWidth);
    exclude=[];
    for j=1:length(locs)-1
       if min(meanRates(i,locs(j):locs(j+1))) > ((pks(j)+pks(j+1))./2) * threshold
           % exclude fields without a 90 % decrease in rate between peaks
           if pks(j) > pks(j+1)
               exclude = [exclude;j+1];
           elseif pks(j) < pks(j+1)
               exclude = [exclude;j];
           end
       end
    end   
    % remove field peaks with a standard dev higher than the mean
    % (unreliable fields)
    exclude = [exclude; find( meanRates(i,locs) <  stdRates(i,locs))'];
    
    pks(exclude) = [];
    locs(exclude)=[];
    fieldCount = 1;
    for j=1:length(locs)
        
        Map_Field = meanRates(i,:) > pks(j) * threshold;
        
        start = locs(j);
        stop = locs(j);
        while Map_Field(start) == 1  && start > 1
            start = start-1;
        end
        while Map_Field(stop) == 1  && stop < length(Map_Field) -1
            stop = stop+1;
        end
        if stop - start > minFieldWidth && stop - start < maxFieldWidth
            fields{i}{fieldCount}.start = start;
            fields{i}{fieldCount}.stop = stop;
            fields{i}{fieldCount}.width = stop - start;
            fields{i}{fieldCount}.peakFR = pks(j);
            fields{i}{fieldCount}.peakLoc = locs(j);
            fields{i}{fieldCount}.spkCount = locs(j);
            com = start; % calculate center of mass for field
            fields{i}{fieldCount}.COM = fields{i}{fieldCount}.peakLoc;
            while sum(meanRates(i,start:stop)) - sum(meanRates(i,start:com)) > sum(meanRates(i,start:stop))./2
                fields{i}{fieldCount}.COM = com;
                com = com + 1;
            end
            fieldCount = fieldCount + 1;
        end
    end
end


warning on



end

