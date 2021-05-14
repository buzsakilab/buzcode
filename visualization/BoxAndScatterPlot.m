function [ datamean,datastd ] = BoxAndScatterPlot( data,varargin )
%BoxAndScatterPlot(data) makes a box plot with data points.
%
%INPUTS
%   data    {ngroups} cell array - each cell has a set of data points
%   (options)
%       'colors'    [ngroupx x 3] matrix of colors for each groups
%       'labels'    {ngroups} cell array of group names
%
%
%DLeventein 2018
%%
if ~iscell(data)
    data = {data};
end
numdatagroups = length(data);

colordefault = zeros(numdatagroups,3);
labelsdefault = 1:numdatagroups;

p = inputParser;
addParameter(p,'colors',colordefault)
addParameter(p,'labels',labelsdefault)
addParameter(p,'groupnumbers',labelsdefault)
addParameter(p,'withpoints',true)
parse(p,varargin{:})
colors = p.Results.colors;
labels = p.Results.labels;
groupnumbers = p.Results.groupnumbers;
withpoints = p.Results.withpoints;

%% Reformat Data, calculate stats
datagrouped = [];
datagroups = [];
for gg = 1:numdatagroups
    if isrow(data{gg})
        data{gg} = data{gg}';
    end
    if isempty(data{gg})
        data{gg} = nan;
    end
    datagrouped = [datagrouped ; data{gg}];
    datagroups = [datagroups ; gg.*ones(size(data{gg}))];
    
    datamean(gg) = mean(data{gg});
    datastd(gg) = std(data{gg});
end
%%

boxplot(datagrouped,datagroups,'symbol','','colors',colors,'labels',labels,'positions',groupnumbers)
hold on
if withpoints
for gg = 1:numdatagroups
    plot(groupnumbers(gg).*ones(size(data{gg}))+0.025.*randn(size(data{gg})),data{gg},'.','color',colors(gg,:))
end
end

end

