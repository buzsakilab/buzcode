function [lfpwhiten] = bz_whitenLFP(lfp,varargin)

% [y, ARmodel] = bz_whitenLFP(lfp,varargin)
% Whiten an lfp signal using an autoregresive model.
% Requires the external package 'arfit'

% INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                   -lfp can also be a [t x 1] timeseries signal. in which
%                   case you need to input 'samplingRate'
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%       window     if window specified will recompute the AR model in each 
%                   window of that size. Default [1 length recording].
%       ARorder    model order. Default 2. 
%       commonAR   if true then will use model from first channel for all,
%                  if false, will calculate one model per channel. Default true.
%       ARmodel    if ARmodel is provided - use it, not compute fromthe data
%                   output optionaly the ARmodel for use on the other data 
%                  to be on the same scale
%           
%    =========================================================================
%
% OUTPUT
%    lfpwhiten            a buzcode structure with fields lfpwhiten.data,
%                                                   lfpwhiten.timestamps
%                                                   lfpwhiten.samplingRate
%                                                   lfpwhiten.params

% AntonioFR, 5/20

%% Parse the inputs

%Parameters
parms = inputParser;
addParameter(parms,'window',[],@isnumeric);
addParameter(parms,'commonAR',true,@islogical);
addParameter(parms,'ARorder',2,@isnumeric);

parse(parms,varargin{:})
window = parms.Results.window;
commonAR = parms.Results.commonAR;
ARorder = parms.Results.ARorder;


%lfp input
if isstruct(lfp)
    samplingRate = lfp.samplingRate;
elseif isempty(lfp)
    wavespec = lfp;
    return
elseif isnumeric(lfp)
    data_temp = lfp;
    clear lfp
    lfp.data = data_temp;
    lfp.timestamps = [1:length(lfp.data)]'./samplingRate;
end

Lin =  lfp.data;
Lout = zeros(size(Lin,1),size(Lin,2));

%%
ARmodel = [];
if isempty(window)
    seg = [1 size(Lin,1)];
else
    nwin = floor(size(Lin,1)/window)+1;
    seg = repmat([1 window],nwin,1)+repmat([0:nwin-1]'*window,1,2);
    if nwin*window > size(Lin,1)
        seg(end,2) = size(Lin,1);
    end   
end

for w=1:size(seg,1)
    if ~isempty(ARmodel) 
        A = ARmodel;
        for i=1:nCh
            Lout(seg(w,1):seg(w,2),i) = filt0(A, Lin(seg(w,1):seg(w,2),i));
        end
    else
        if commonAR  
            for i=1:size(Lin,2)
                if  w==1 && i==1
                    [k,Atmp] = arfit(Lin(seg(w,1):seg(w,2),i),ARorder,ARorder);
                    A = [1 -Atmp];
                    ARmodel = A;
                end
                Lout(seg(w,1):seg(w,2),i) = filt0(Lin(seg(w,1):seg(w,2),i),A);
            end
        else
            for i=1:size(Lin,2)
                switch artype
                    case 1
                        [k,Atmp] = arfit(Lin(seg(w,1):seg(w,2),i),ARorder,ARorder);
                        A =[1 -Atmp];
                    case 2
                        A = arburg(Lin(seg(w,1):seg(w,2),i),ARorder);
                end
                Lout(seg(w,1):seg(w,2),i) = filt0(Lin(seg(w,1):seg(w,2),i), A);
            end
        end
    end
end

%% output 

lfpwhiten.data = Lout;
lfpwhiten.timestamps = lfp.timestamps;
lfpwhiten.samplingRate = lfp.samplingRate;
lfpwhiten.params.whitening = true;
lfpwhiten.params.ARorder = ARorder;
lfpwhiten.params.commomAR = commonAR;
lfpwhiten.params.ARmodel = A;

end


%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = filt0(x,W)

if all(size(W)>=2), error('window must be a vector'), end
if numel(x)==max(size(x)), x=x(:); end

C = length(W);
if C > size( x, 1 )
    Y = NaN * ones( size( x ) );
    return
end
D = ceil(C/2) - 1;
Y = filter(W,1,[flipud(x(1:C,:)); x; flipud(x(end-C+1:end,:))]);
Y = Y(1+C+D:end-C+D,:);

end
