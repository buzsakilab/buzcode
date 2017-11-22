function mono_res = bz_MonoSynConvClick (spikeIDs,spiketimes,varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  INPUTS
%%%
%%%  spikeIDs = Nx3 matrix. N = # spikes.
%%%      col1 = shank ID, col2 = cluID on shanke, col3 = unit ID
%%%
%%%  spikestimes = Nx1 matrix aligned with spikeIDs where each row is a
%%%      time stamp in ms.
%%%
%%%  OPTIONAL INPUTS:
%%%
%%%  binSize = timebin to compute CCG (in seconds)
%%%
%%%  duration = window to compute CCG (in seconds)
%%%
%%%  epoch = [start end] (in seconds)
%%%
%%%  cells = N x 2 matrix -  [sh celID] to include (NOTE indexing will be
%%%          done on full spikeIDlist
%%%
%%%  conv_w = # of time bins for computing the CI for the CCG.
%%%
%%%  alpha = type I error for significance testing using convolution
%%%     technique. Stark et al, 2009
%%%
%%%  calls: CCG, InInterval,FindInInterval (from FMA toolbox)
%%%         tight_subplot, mtit (from matlabcentral)
%%%         bz_cch_conv, bz_PlotMonoSyn
%%%
%%%  OUTPUT
%%%  mono_res.alpha = p-value
%%%  mono_res.ccgR = 3D CCG (time x ref x target;
%%%  mono_res.sig_con = list of significant CCG;
%%%  mono_res.Pred = predicted Poisson rate;
%%%  mono_res.Bounds = conf. intervals of Poisson rate;
%%%  mono_res.conv_w = convolution windows (ms)
%%%  mono_res.completeIndex = cell ID index;
%%%  mono_res.binSize = binSize;
%%%  mono_res.duration = duration;
%%%  mono_res.manualEdit = visual confirmation of connections
%%%  mono_res.Pcausal = probability of getting more excess in the causal than anticausal direction;
%%%  mono_res.FalsePositive = FalsePositive rate from English et al., 2017;
%%%  mono_res.TruePositive = TruePositive rate from English et al., 2017;

%%%  EXAMPLE:
%%%
%%%  mono_res = bz_MonoSynConvClick (spikesIDs,spiketimes,'binsize',.0005,'duration',.2, ...
%%%  'alpha',.05,'conv_w',20,'cells',[1 2;1 3;4 5;8 8],'epoch',[10 3000; 4000 5000]);
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%


%get experimentally validated probabilities


fil = which('bz_MonoSynConvClick');
sl = regexp(fil,'/');
fil = fil(1:sl(end));
  foundMat = false;
if exist([fil 'ProbSynMat.mat'],'file')==2
    v = load([fil 'ProbSynMat.mat']);
    foundMat = true;
else
    warning('You do not have the ProbSynMat matrix describing the likelihood of experimentally validated connectivty given excess syncrony')
end


%parse inputs

%set defaults
binSize = .0004; %.4ms
duration = .2; %200ms
epoch = [0 inf]; %whole session
cells = unique(spikeIDs(:,1:2),'rows');
nCel = size(cells,1);
conv_w = .010/binSize;  % 10ms window
alpha = 0.001; %high frequency cut off, must be .001 for causal p-value matrix
plotit = true;
sorted = false;

if length(varargin) ==1 && iscell(varargin{1})
    varargin = varargin{1};
end


% Parse options
for i = 1:2:length(varargin),
    if ~isa(varargin{i},'char'),
        error(['Parameter ' num2str(i+3) ' is not a property ']);
    end
    switch(lower(varargin{i})),
        case 'duration',
            duration = varargin{i+1};
            if ~isa(duration,'numeric') | length(duration) ~= 1 | duration < 0,
                error('Incorrect value for property ''duration''');
            end
        case 'binsize',
            binSize = varargin{i+1};
            if ~isa(binSize,'numeric') | length(binSize) ~= 1 | binSize <= 0,
                error('Incorrect value for property ''binsize'' ');
            end
        case 'epoch',
            epoch = varargin{i+1};
            if ~isa(epoch,'numeric') | size(epoch,2) ~= 2,
                error('Incorrect value for property ''epoch'' ');
            end
            
        case 'cells',
            cells = varargin{i+1};
          
        case 'conv_w',
            conv_w = varargin{i+1};
           
        case 'alpha',
            alpha = varargin{i+1};
            if ~isa(alpha,'numeric'),
                error('Incorrect value for property ''alpha''');
            end
            
        case 'plot',
            plotit = varargin{i+1};
            if ~islogical(plotit),
                error('Incorrect value for property ''plot''');
            end
            
        case 'sorted'
            sorted = varargin{i+1};
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help LoadBinary'' for details).']);
    end
end

if ~sorted
    %sort by spike times
    [spiketimes,b] = sort(spiketimes);
    spikeIDs = spikeIDs(b,:);
end

%restrict by cells and epochs

[status] = InIntervals(spiketimes,epoch);
allID = unique(spikeIDs(:,3));
kp = ismember(spikeIDs(:,1:2),cells,'rows') & status;
spikeIDs = spikeIDs(kp,:);
spiketimes = spiketimes(kp);

nBins=duration/binSize+1;
[IDindex,tet_idx,ID_idx] = unique(spikeIDs(:,3));	% list of IDs

completeIndex = spikeIDs(tet_idx,:);


% Create CCGs (including autoCG) for all cells
[ccgR1,tR] = CCG(spiketimes,spikeIDs(:,3),'binSize',binSize,'duration',duration);

ccgR = nan(size(ccgR1,1),nCel,nCel);
ccgR(:,1:size(ccgR1,2),1:size(ccgR1,2)) = ccgR1;


% get  CI for each CCG
Pval=nan(length(tR),nCel,nCel);
Pred=zeros(length(tR),nCel,nCel);
Bounds=zeros(size(ccgR,1),nCel,nCel);
sig_con = [];
sig_con1 = [];

TruePositive = nan(nCel,nCel);
FalsePositive = nan(nCel,nCel);
Pcausal = nan(nCel,nCel);
for refcellID=1:max(IDindex)
    for cell2ID= refcellID+1:max(IDindex)
        
        cch=ccgR(:,refcellID,cell2ID);			% extract corresponding cross-correlation histogram vector
        
        refcellshank=completeIndex(completeIndex(:,3)==refcellID);
        cell2shank=completeIndex(completeIndex(:,3)==cell2ID);
        if refcellshank==cell2shank
            
            
            % central 1.6 ms on same-shank = NaN due to limitations of
            % extracting overlapping spikes
            
            centerbins = ceil(length(cch)/2);
            sameshankcch=cch;
            sameshankcch(centerbins)=[];
            
            [pvals,pred,qvals]=bz_cch_conv(sameshankcch,conv_w);
            pred=[pred(1:(centerbins(1)-1));nan(length(centerbins),1);pred(centerbins(end)-length(centerbins)+1:end)];
            
            pvals=[pvals(1:(centerbins(1)-1));nan(length(centerbins),1);pvals(centerbins(end)-length(centerbins)+1:end)];
        else
            % calculate predictions using Eran's bz_cch_conv
            
            [pvals,pred,qvals]=bz_cch_conv(cch,conv_w);
        end
        
        
        
        % Store predicted values and pvalues for subsequent plotting
        Pred(:,refcellID,cell2ID)=pred;
        Pval(:,refcellID,cell2ID)=pvals(:);
        Pred(:,cell2ID,refcellID)=flipud(pred(:));
        Pval(:,cell2ID,refcellID)=flipud(pvals(:));
        
        % Calculate upper and lower limits with bonferonni correction
        % monosynaptic connection will be +/- 4 ms
        
        nBonf = round(.004/binSize)*2;
        
        hiBound=poissinv(1-alpha/nBonf,pred);
        loBound=poissinv(alpha/nBonf, pred);
        Bounds(:,refcellID,cell2ID,1)=hiBound;
        Bounds(:,refcellID,cell2ID,2)=loBound;
        
        Bounds(:,cell2ID,refcellID,1)=flipud(hiBound(:));
        Bounds(:,cell2ID,refcellID,2)=flipud(loBound(:));
        
        % sig = cch>hiBound | cch < loBound;
        sig = cch>hiBound;
        
        % Find if significant periods falls in monosynaptic window +/- 4ms
        prebins = round(length(cch)/2 - .0032/binSize):round(length(cch)/2);
        postbins = round(length(cch)/2 + .0008/binSize):round(length(cch)/2 + .004/binSize);
        cchud  = flipud(cch);
        sigud  = flipud(sig);
        sigpost=max(cch(postbins))>poissinv(1-alpha,max(cch(prebins)));
        sigpre=max(cchud(postbins))>poissinv(1-alpha,max(cchud(prebins)));
        
        
        %define likelihood of being a connection
        pvals_causal = 1 - poisscdf( max(cch(postbins)) - 1, max(cch(prebins) )) - poisspdf( max(cch(postbins)), max(cch(prebins)  )) * 0.5;
        pvals_causalud = 1 - poisscdf( max(cchud(postbins)) - 1, max(cchud(prebins) )) - poisspdf( max(cchud(postbins)), max(cchud(prebins)  )) * 0.5;
        
        %can go negative for very small p-val - beyond comp. sig. dig
        
        if pvals_causalud<0
            pvals_causalud = 0;
        end
        
        if pvals_causal<0
            pvals_causal = 0;
        end
        
        
        Pcausal(refcellID,cell2ID) = pvals_causal;
        Pcausal(cell2ID,refcellID) = pvals_causalud;
        if foundMat
            if any(Pval(postbins,cell2ID,refcellID)<.001)
                
                FP =  v.ProbSyn.FalsePositive((histc(pvals_causalud,v.ProbSyn.thres))>0);
                TP =  v.ProbSyn.TruePositive((histc(pvals_causalud,v.ProbSyn.thres))>0);
                TruePositive(cell2ID,refcellID) = TP;
                FalsePositive(cell2ID,refcellID) = FP;
            end
            
            
            if any(Pval(postbins,refcellID,cell2ID)<.001)
                
                FP =  v.ProbSyn.FalsePositive((histc(pvals_causal,v.ProbSyn.thres))>0);
                TP =  v.ProbSyn.TruePositive((histc(pvals_causal,v.ProbSyn.thres))>0);
                TruePositive(refcellID,cell2ID) = TP;
                FalsePositive(refcellID,cell2ID) = FP;
            end
        end
        
        %check which is bigger
        if (any(sigud(prebins)) && sigpre)
            
            %test if causal is bigger than anti causal
            
            
            sig_con = [sig_con;cell2ID refcellID];
            
        end
        
        if any(sig(postbins)) && sigpost
            
            sig_con = [sig_con;refcellID cell2ID];
        end
        
        
        
    end
    
end

%plot
if plotit
    sig_con = bz_PlotMonoSyn(ccgR,sig_con,Pred,Bounds,completeIndex,binSize,duration);
end
nCel = size(completeIndex,1);
n = histc(spikeIDs(:,3),1:length(allID));
[nn1,nn2] = meshgrid(n);

temp = ccgR - Pred;
prob = temp./permute(repmat(nn2,1,1,size(ccgR,1)),[3 1 2]);



%save outputs
mono_res.ccgR = ccgR;
mono_res.Pval = Pval;
mono_res.prob = prob;
mono_res.prob_noncor = ccgR./permute(repmat(nn2,1,1,size(ccgR,1)),[3 1 2]);
mono_res.n = n;
mono_res.sig_con = sig_con;
mono_res.Pred = Pred;
mono_res.Bounds = Bounds;
mono_res.completeIndex = completeIndex;
mono_res.binSize = binSize;
mono_res.duration = duration;
%mono_res.spiketimes = spiketimes;
%mono_res.spikeIDs = spikeIDs;
mono_res.ManuelEdit = plotit;
mono_res.conv_w = conv_w;
mono_res.Pcausal = Pcausal;


if foundMat
mono_res.FalsePositive = FalsePositive;
mono_res.TruePositive = TruePositive;

end

end
