function [ CONDXY ] = ConditionalHist( X,Y,varargin)
%[CONDXY] = ConditionalHist(X,Y) Calculates the conditional probabilty 
%of Y given X
%
%INPUT
%   X
%   Y
% (options)
%   'numXbins'
%   'numYbins'
%   'Xbounds'
%   'Ybounds'
%   'minX'
%   
%
%OUTPUT
%   CONDXY
%       .pYX    a [numXbins x numYbins] matrix in which the [x,y]th element is P(y|x) 
%       .XYhist joint histogram of X and Y
%       .XYprop joint probability of X and Y
%       .Xhist  histogram of X
%       .Xbins  X bins
%       .Ybins  Y bins
%
%DLevenstein 2019
%%
p = inputParser;
addParameter(p,'numXbins',50)
addParameter(p,'numYbins',50)
addParameter(p,'Xbounds',[])
addParameter(p,'Ybounds',[])
addParameter(p,'minX',25)
parse(p,varargin{:})
numXbins = p.Results.numXbins;
numYbins = p.Results.numYbins;
Xbounds = p.Results.Xbounds;
Ybounds = p.Results.Ybounds;
minX = p.Results.minX;


%% For cell input

if iscell(Y) && iscell(X)
    CONDXY = cellfun(@(x,y) ConditionalHist(x,y,varargin{:}),...
        X,Y,'UniformOutput',false);
    CONDXY = bz_CollapseStruct([CONDXY{:}],3);
    CONDXY.Xbins = CONDXY.Xbins(1,:,1);
    CONDXY.Ybins = CONDXY.Ybins(1,:,1);
    return
end

%% For multiple columns in X - conditonal probabilty of each

if size(Y,2)>1 && size(X,2)==1
    for yy = 1:size(Y,2)
        CONDXY(yy) = ConditionalHist(X,Y(:,yy),varargin{:});
    end
    CONDXY = bz_CollapseStruct(CONDXY,3);
    return
end



%%
if isempty(Xbounds)
    Xbounds(1) = min(X); Xbounds(2) = max(X);
end
if isempty(Ybounds)
    Ybounds(1) = min(Y(~isinf(Y))); Ybounds(2) = max(Y(~isinf(Y)));
end

Xedges = linspace(Xbounds(1),Xbounds(2),numXbins+1);
Xbins = Xedges(1:end-1)+ 0.5.*diff(Xedges([1 2]));
Xedges(1) = -inf;Xedges(end) = inf;

Yedges = linspace(Ybounds(1),Ybounds(2),numYbins+1);
Ybins = Yedges(1:end-1)+ 0.5.*diff(Yedges([1 2]));
Yedges(1) = -inf;Yedges(end) = inf;

%First calculate the marginal probability of X
[Xhist,~,XbinID] = histcounts(X,Xedges);
Xhist4norm = Xhist;Xhist4norm(Xhist4norm<=minX) = nan;

%Then calculate the joint probabilty of X and Y
if length(Ybins) ==1
    XYhist = Xhist;
else
    [XYhist] = hist3([X,Y],{Xbins,Ybins});
end

% Conditional probability of Y given X
pYX = bsxfun(@(x,y) x./y,XYhist,Xhist4norm');

%Mean Y given X
for xx = 1:length(Xbins)
    meanYX(xx) = nanmean(Y(XbinID==xx));
    meanYX(isnan(Xhist4norm)) = nan;
end

CONDXY.pYX = pYX;
CONDXY.XYhist = XYhist;
CONDXY.XYprob = XYhist./nansum(XYhist(:));
CONDXY.meanYX = meanYX;
CONDXY.Xhist = Xhist;
CONDXY.Xbins = Xbins;
CONDXY.Ybins = Ybins;


%%
% figure
% imagesc(CONDXY.Xbins,CONDXY.Ybins,CONDXY.pYX')
% axis xy
end

