%
%  function [ccgR,tR,GSPExc,GSPInh,ccgjMtx,ccgjstats]=CCG_jitter(spiket,spikeind,clu1,clu2,dataSampFreq,BinSize,HalfBins,jscale,njitter,alpha,PresenBoolen)
%
%  ---------------------------------------------------
%  INPUTS
%  spiket, spikeind  : res file, clu file
%  clu1, clu2        : the cell numbers
%  dataSampFreq      : original sampling rate of the data in Hz (usually 20000) 
%  Bizsize, HalfBins : --> see 'CCG.m'
%  jscale            : jittering scale, unit is 'ms'
%  njitter           : # of jittering
%  alpha             : significant level... measured in occurence frequency
%  ---------------------------------------------------
%  OUTPUTS
%  ccgR                   : ccg of real data   <-- [ccg,t]=CCG(...);
%  tR                     : timepoints of real data
%  GSPExc                 : Global Significant Period of Mono Excitation.
%  GSPInh                 : Global Significant Period of Mono Inhibition.
%  ccgjMtx                : Jitter matrix
%  ccgjstats.pointwiseMax : Signficance cutoff from pointwise max jitter
%  ccgjstats.pointwiseMin : Signficance cutoff from pointwise min jitter
%  ccgjstats.globalMax    : Signficance cutoff from global max jitter
%  ccgjstats.globalMin    : Signficance cutoff from global min jitter
%  ---------------------------------------------------
%  Example Inputs: CCG_jitter(spiket,spikeind,clu1,clu2,20000,20,40,2,500,0.01,1);
%             --- 20000Hz sample rate, binsize is 1ms (20samples), 40 bins,
%             jitter time scale is 2ms, 500 times jittering, p<0.01, plot
%             output in figure
%
%  # Sampling rate was initially fixed to 20khz, should now be generalized
%      - BW
%
%  Coded by  Shigeyoshi Fujisawa
%   based on Asohoan Amarasingham's resampling methods
%  Very slight modifications/generalizations by Brendon Watson

function [ccgR,tR,GSPExc,GSPInh,ccgjMtx,ccgjstats]=CCG_jitter(spiket,spikeind,clu1,clu2,dataSampFreq,BinSize,HalfBins,jscale,njitter,alpha,PresenBoolen)

if nargin <11
   PresenBoolen=1; % make figure or not
end

   res1 = spiket(spikeind==clu1);
   res2 = spiket(spikeind==clu2);

   [ccgR, tR] = CCG([res1;res2],[ones(size(res1));2*ones(size(res2))], BinSize, HalfBins, dataSampFreq,[1,2],'count');

   %%%%%%  CCG for jittering data
   for i=1:njitter
%        res2_jitter = res2 + 2*(20*jscale)*rand(size(res2))-1*20*jscale; 
       jitterms = round(dataSampFreq/1000*jscale);
       res2_jitter = res2 + 2*(jitterms)*rand(size(res2))-1*jitterms; 
       [ccg, tJ] = CCG([res1;res2_jitter],[ones(size(res1));2*ones(size(res2))], BinSize, HalfBins, dataSampFreq,[1,2],'count');
       ccgj(:,i)=ccg(:,1,2);
       ccgjmax(i)=max(ccgj(:,i));
       ccgjmin(i)=min(ccgj(:,i));
   end

   %%%%%%  Compute the pointwise line
   signifpoint = njitter*alpha;
   for i=1:length(tJ)
       sortjitterDescend  = sort(ccgj(i,:),'descend');
       sortjitterAscend   = sort(ccgj(i,:),'ascend');
       ccgjptMax(i) = sortjitterDescend(signifpoint);
       ccgjptMin(i) = sortjitterAscend(signifpoint);
   end

   %%%%%%  Compute the global line
   sortgbDescend   = sort(ccgjmax,'descend');
   sortgbAscend    = sort(ccgjmin,'ascend');
   ccgjgbMax  = sortgbDescend(signifpoint)*ones(size(tJ));
   ccgjgbMin  = sortgbAscend(signifpoint)*ones(size(tJ));

   ccgjm  = mean(ccgj,2);


%%%%%%%%%%%%%% Significant Period

if nargout>2

   findExc = find((ccgR(:,1,2)>ccgjgbMax')&(ccgR(:,1,2)>0));
   findInh = find((ccgR(:,1,2)<ccgjgbMin')&(ccgjgbMin'>0));

   GSPExc = zeros(size(tR));  % Global Significant Period of Mono Excitation
   GSPInh = zeros(size(tR));  % Global Significant Period of Mono Inhibition

   GSPExc(findExc) = 1;
   GSPInh(findInh) = 1;

end


%%%%%%%%%%%%%%% Prep for output
if nargout>4
   ccgjMtx=ccgj;
end

ccgjstats.pointwiseMax = ccgjptMax;
ccgjstats.pointwiseMin = ccgjptMin;
ccgjstats.globalMax = ccgjgbMax;
ccgjstats.globalMin = ccgjgbMin;



%%%%%%%%%%%%%%% Presentation, if chosen

if PresenBoolen==1
      %plot(tR,ccgR(:,1,2),'color','k')
      bar(tR,ccgR(:,1,2),'facecolor','k','edgecolor','k')
      line(tR,ccgjm,'linestyle','--','color','b')
      line(tR,ccgjptMax,'linestyle','--','color','r')
      line(tR,ccgjgbMax,'linestyle','--','color','m')
      line(tR,ccgjptMin,'linestyle','--','color','r')
      line(tR,ccgjgbMin,'linestyle','--','color','m')
      set(gca,'XLim',[min(tR),max(tR)])
end

if 0
   % if you want to plot y as 'hz'
   [ccgR_hz, tR] = CCG([res1;res2],[ones(size(res1));2*ones(size(res2))], BinSize, HalfBins, dataSampFreq,[1,2],'hz');
   modif = sum(ccgR_hz(:,1,2))/sum(ccgR(:,1,2));
   modif = 1 % 'count'
   %
   bar(tR,modif*ccgR(:,1,2),'facecolor','k','edgecolor','k')
   line(tR,modif*ccgjm,'linestyle','--','color','b')
   line(tR,modif*ccgjptMax,'linestyle','--','color','r')
   line(tR,modif*ccgjgbMax,'linestyle','--','color','m')
   line(tR,modif*ccgjptMin,'linestyle','--','color','r')
   line(tR,modif*ccgjgbMin,'linestyle','--','color','m')
   set(gca,'XLim',[min(tR),max(tR)])
end