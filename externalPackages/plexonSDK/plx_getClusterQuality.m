function [LRatio IsolDist ID] = plx_getClusterQuality(varargin)
% USAGE
% 
% [LRatio IsolDist ID] = getClusterQuality(plxfile)
%
% INPUTS
%       
%  plx - filename of a *.plx file in the current working directory (should
%        already be manually cluster cut)
%
% OUTPUTS
%
%  LRatio - L_Ratio.m
%  IsolDist - IsolationDistance.m
%  ID - 
%
%
% pulls the spike waveforms out of a .plx file and returns cluster quality
% metrics

if nargin < 1
    plx = dir('*plx'); 
    for i=1:length(plx)
        if ~isempty(findstr(lower(plx(i).name),'cut'))
            plx = plx(i).name;
        end
    end
else
    plx = varargin{1};    
end

if isempty(plx)
    LRatio=nan;
    IsolDist=nan;
    ID=nan;
    return 
end
try
    [OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreTresh, ...
        SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, ...
        DateTime] = plx_information(plx);


    for i=1:Trodalness:48
        for j=1:50  % ..if you're recording > 50 units per *trode, you may want to do some merging...
        [n, npw, ts, wave{i,j}] = plx_waves_v(plx,i,j-1); % -1 includes 'unsorted spikes'
        [n, npw, ts, wave{i+1,j}] = plx_waves_v(plx,i+1,j-1);
        end
    end

    c=1;
    for j=1:Trodalness:48
        w=[];id=[];
        for i=1:50
            if wave{j,i}~=-1
                if size(wave{j,i},1) == size(wave{j+1,i},1)
                w =[w; wave{j,i}, wave{j+1,i}];
                id = [id;i*ones(size(wave{j,i},1),1)];
                end
            end
        end
        for i=1:50
            if ~isempty(find(id==i))
           IsolDist(c) = IsolationDistance(w, find(id==i));
           [L, LRatio(c), df]  = L_Ratio(w, find(id==i));
           if i == 1
               ID(c) = 0;
           else
               ID(c) = 1;
           end
           c=1+c;
            end
        end
        clear w id
    end
end


