function [normdata,intmean,intstd] = NormToInt(data,normtype,int,sf,varargin)
%NormToInt(data,normtype, int,sf) normalizes the data to a subset of time 
%(int). Has options for multiple normalization types.
%
%INPUTS
%   data    [Nt x Ndim] 
%   normtype 'mean','median' 'Z', 'max', 'percentile', 'modZ' for modified Z-score
%   int     [Nints x 2] reference interval onsets and offset times   
%           to normalize the data with respect to (optional)
%   sf      (optional) sampling frequency of the data. default 1
%
%   (options)
%   'moving'    window for moving normatilization
%   'SHOWFIG'   show a figue with distribution of the data
%Note: modified Z score is median-based and robust to outliers
%see https://ibm.co/2qi4Vy5
%
%DLevenstein Summer 2016
%TO DO: Improve input parsing, compadible with buzcode structures
%%
p = inputParser;
addParameter(p,'moving',0,@isnumeric);
addParameter(p,'SHOWFIG',false);
parse(p,varargin{:})
movingwin = p.Results.moving;
SHOWFIG = p.Results.SHOWFIG;
if movingwin>0
    MOVING = true;
else
    MOVING = false;
end
%%
if ~exist('int','var') || isempty(int)
    int = [1 size(data,1)];
end

if isa(int,'intervalSet')
    int = [Start(int,'s'), End(int,'s')];
end


numints = length(int(:,1));
if ~exist('sf','var')|| isempty(sf)
    sf = 1;
end


int = round(int*sf); %Convert intervals from seconds to indices 
int(int==0)=1; %Turn time=0 to the first indexed datapoint
int(isinf(int))=size(data,1);

%Make vector that is nans at times not in the intervals
int_data = nan(size(data));
for ii = 1:numints
    int_data(int(ii,1):int(ii,2),:) = data(int(ii,1):int(ii,2),:);
end 

switch MOVING
    case false        
        switch normtype
            case {'Z','mean'}
                intmean = nanmean(int_data,1);
                intstd = nanstd(int_data,0,1);
                intmeanlong = repmat(intmean,length(data(:,1)),1);
                intstdlong = repmat(intstd,length(data(:,1)),1);
            case {'modZ','median'}
                intmedian = nanmedian(int_data);
                intMAD = mad(int_data,[],1);
                intmedianlong = repmat(intmedian,length(data(:,1)),1);
                intMADlong = repmat(intMAD,length(data(:,1)),1);
        end
    case true
        switch normtype
            case {'Z','mean'}
                intmean = movmean(int_data,movingwin.*sf,1,'omitnan');
                intstd = movstd(int_data,movingwin.*sf,'omitnan');
                %Because moving mean will have overhangs out of int
                intmean(isnan(int_data))=nan; intstd(isnan(int_data))=nan; 
                intmeanlong = intmean; intstdlong = intstd;

            case {'modZ','median'}
                intmedian = movmedian(int_data,movingwin.*sf,1,'omitnan');
                intMAD = movmad(int_data,movingwin.*sf,'omitnan');
                intmedian(isnan(int_data))=nan; intMAD(isnan(int_data))=nan; 
                intmedianlong = intmedian; intMADlong = intMAD;

        end
end        

switch normtype
    case 'Z'
        normdata = (data-intmeanlong)./intstdlong;
    case 'mean'
        normdata = data./intmeanlong;
    case 'median'
        normdata = data./intmedianlong;
    case 'max'
        colmax = max(int_data,[],1);
        normdata = bsxfun(@(X,Y) X./Y,data,colmax);
    case 'percentile'
        sortdata = unique(sort(int_data(~isnan(int_data))));
        percentiles = linspace(0,1,length(sortdata));
        normdata = interp1(sortdata,percentiles,data,'nearest');
    case 'modZ'
        normdata = 0.6745*(data-intmedianlong)./intMADlong;
        intmean = intmedian;
        intstd = intMAD;
    otherwise
        display('incorrect normtype')
end


%%
if SHOWFIG
    figure
        subplot(2,1,1)
            plot(data,normdata,'.')
        subplot(2,1,2)
            hist(int_data)
end


end

