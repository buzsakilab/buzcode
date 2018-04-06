function [thresh,cross,bihist,diptest] = bz_BimodalThresh(bimodaldata,varargin)
%[thresh,cross,bihist] = BimodalThresh(bimodaldata) takes bimodal time
%series data, calculates the threshold between the modes (i.e. UP vs DOWN states),
%and returns the crossing times (i.e. UP/DOWN onset/offset times)
%
%INPUTS
%   bimodaldata     vector of bimodal data
%
%       (optional)  for optional inputs, use: ('inputoption', inputvalue)
%       'maxthresh' sets a maximum threshold
%       'Schmidt'   Schmidt trigger uses halfway points between trough and 
%                   lower/upper peaks for DOWN/UP state thresholds
%       'maxhistbins' Maximum number of hist bins to try before giving up
%       'startbins'  minimum number of hist bins (initial hist)
%       'setthresh'     set your own threshold
%
%OUTPUS
%   thresh          threshold
%   cross
%       .upints     intervals of UP states [onsets offsets]
%       .downints   intervals of DOWN states [onsets offsets]
%   bihist
%       .bins       bin centers
%       .hist       counts
%   diptest
%       .dip             hartigans dip test statistic
%       .p              p value for bimodal distribution
%
%DLevenstein Spring 2016
%-Updated Spring 2017 to include Schmidt trigger, maxthresh
%-TO DO: add hartigans dip test to determine if bimodal... also to find the
%best threshold?  Or ISO-CUT from Simons foundation?
%%
%Optional input defaults
maxthresh = inf;
Schmidt = false;
maxhistbins = 25;
SETTHRESH = false;
starthistbins = 10;

for i = 1:length(varargin)
	switch(lower(varargin{i}))
        case 'maxthresh'
            maxthresh = varargin{i+1};
        case 'maxhistbins'
            maxhistbins = varargin{i+1};
        case 'startbins'
            starthistbins = varargin{i+1};
        case 'schmidt'
            Schmidt = true;
        case 'setthresh'
            thresh = varargin{i+1};
            SETTHRESH = true;
    end
end

%Run hartigansdiptest for bimodality
nboot = 500; %number of times for bootstrapped significance
[diptest.dip, diptest.p_value]=hartigansdipsigniftest(bimodaldata,nboot);

if diptest.p_value>0.05  %add option to bipass this
    display('Dip test says: not bimodal')
    cross.upints = []; cross.downints = []; thresh=nan; 
    [bihist.hist,bihist.bins]= hist(bimodaldata,starthistbins);
    return
end

%Remove data over the threshold... this is klugey
overmax = bimodaldata>=maxthresh;
bimodaldata(overmax) = nan;

numpeaks = 1;
numbins = starthistbins; 
while numpeaks ~=2
    [bihist.hist,bihist.bins]= hist(bimodaldata,numbins);

    [PKS,LOCS] = findpeaks([0 bihist.hist 0],'NPeaks',2);
    LOCS = sort(LOCS)-1;
    numbins = numbins+1;
    numpeaks = length(LOCS);
    if numbins==maxhistbins && ~SETTHRESH
        display('Unable to find trough')
    	cross.upints = []; cross.downints = []; thresh=nan; return
    end
end

betweenpeaks = bihist.bins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks(-bihist.hist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

if ~SETTHRESH
thresh = betweenpeaks(diploc);
end

    




%%

%The Schmidt trigger uses thresholds halfway between trough and peak
%on each end
if Schmidt
    threshUP = thresh + 0.5.*(betweenpeaks(end)-thresh);
    threshDOWN = thresh + 0.5.*(betweenpeaks(1)-thresh);
    
    overUP = bimodaldata>threshUP;
    overDOWN = bimodaldata>threshDOWN;
    
    crossup = find(diff(overUP)==1);
    crossdown = find(diff(overDOWN)==-1);
    
    %Delete incomplete (repeat) crossings
    allcrossings = [crossup ones(size(crossup)) ;...
        crossdown zeros(size(crossdown))];
    [~,sortorder] = sort(allcrossings(:,1));
    allcrossings = allcrossings(sortorder,:);
    updownswitch = diff(allcrossings(:,2));
    samestate = find(updownswitch==0)+1;
    allcrossings(samestate,:)=[];
    
    crossup = allcrossings(allcrossings(:,2)==1,1);
    crossdown = allcrossings(allcrossings(:,2)==0,1);
    
else
    overind = bimodaldata>thresh;
    crossup = find(diff(overind)==1);
    crossdown = find(diff(overind)==-1);
end

%If there's only one crossing...
if isempty(crossup) || isempty(crossdown)
    cross.upints = []; cross.downints = []; return
end

%%
upforup = crossup;
upfordown = crossup;
downforup = crossdown;
downfordown = crossdown;

if crossup(1) < crossdown(1)
    upfordown(1) = [];
end
if crossdown(end) > crossup(end)
    downfordown(end) = [];
end
if crossdown(1) < crossup(1)
    downforup(1) = [];
end
if crossup(end) > crossdown(end)
    upforup(end) = [];
end
    


cross.upints = [upforup downforup];
cross.downints = [downfordown upfordown];


end

