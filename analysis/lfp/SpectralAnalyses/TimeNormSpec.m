function [meaneventspec,t_scaled,freqs,alleventspec,alleventLFP] = TimeNormSpec(LFP,ints,tbins,frange,nfreqs,ncyc,sf,restrictint,intype)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%Inputs:
%   LFP signal
%   ints        [nInts x 2] matrix of interval start and end times
%   tbins       number of time bins
%   
%
%   intype      add option to input already spectrogramed LFP (or
%               distribution or other matrix input signal)
%
%
%Last Updated
%8/15/15
%DLevenstein


%Function Inputs
% intstarts = isse.intstarts ;
% intends = isse.intends;

%Make Defaults
% tbins = 150;
% frange = [1 128];
% nfreqs = 100;   
% ncyc = 10;
% sf = 1250;


goodtime = [Start(restrictint,'s') End(restrictint,'s')];
intstarts = ints(:,1);
intends = ints(:,2);


%%



numints = length(intstarts);

intstarts = intstarts(intstarts>goodtime(1) & intstarts<=goodtime(2));
intends = intends(intends>goodtime(1) & intends<=goodtime(2));

if intends(1)<intstarts(1)
    intends(1) = [];
end
if intends(end)<intstarts(end)
    intstarts(end) = [];
end

numints = length(intstarts);


%% Spectrogram


if sum(size(LFP)==1)==1
    %Downsample LFP - make this the largest possible for 2x max freq
    min_sf = 2*frange(2);
    downsamplefactor = floor(sf/min_sf);

    if downsamplefactor ~= 0
        sf = sf/downsamplefactor;
        LFP = downsample(LFP,downsamplefactor);
    end
    t_LFP = [1:length(LFP)]./sf;

    LFP = LFP(t_LFP>=goodtime(1) & t_LFP<=goodtime(2));
    t_LFP = t_LFP(t_LFP>=goodtime(1) & t_LFP<=goodtime(2));

    display('Calculating Wavelet Spectrogram...')

    [ freqs, t_spec, spec ] = WaveSpec(LFP,frange,nfreqs,ncyc,1/sf,'log');

    rawspec = spec;
    spec = abs(spec);
    zspec = zscore(spec,[],2);
elseif sum(size(LFP)==1)==0
    spec= LFP;
    nfreqs = length(LFP(:,1));
else
    display('Predicted bug detected....')
end

%also look at phase.... which phases are reliable and which have high
%variance?

%convert UPstarts to sample time
UPstarts_si = round(intstarts*sf-goodtime(1)*sf);
UPends_si = round(intends*sf-goodtime(1)*sf);

%%
% [talign,spec_aligned] = AlignToTrigger(zspec,UPstarts_si,W,1/sf_LFP,'spec');
    
    
%% 
display('Scaling Interval Spectrograms...')

numtbins = 3*tbins; %more bins on each side!
t_scaled = linspace(-1,2,numtbins);

alleventspec = zeros(numtbins,nfreqs,numints);
alleventLFP = zeros(numtbins,numints);
newspec = [1 1];    %Placeholder
for u = 1:numints;
    %This is kind of silly/roundabout... interval window with equal window
    %on each size
    ustart = UPstarts_si(u);
    uend = UPends_si(u);
    ulen = uend-ustart+1;
    ustart = ustart-ulen;
    uend = uend+ulen;
    ulen = uend-ustart+1;
    
    %FIX THIS!!!!
    if ustart<0
%         alleventspec(:,:,u) = [];
%         alleventLFP(:,u) = [];
        continue
        
    end
    
    if uend>length(spec(1,:))
%         alleventspec(:,:,u) = [];
%         alleventLFP(:,u) = [];
        continue
    end
    
    uspec = spec(:,ustart:uend)';
    uLFP = LFP(ustart:uend);
    
    
    resamplefact = numtbins/ulen;
    tol = 0.0001;
    while length(newspec(:,1)) ~= numtbins
        [P,Q] = rat(resamplefact,tol);
        if P==0
            tol = tol/10;
            continue
        end
        newspec = resample(uspec,P,Q);
        newLFP = resample(uLFP,P,Q);
        tol = tol/10;
    end
%     bintrate = sf_LFP*P/Q;
%     addedtime = ceil(extrabins/bintrate)*sf_LFP;
%     
%     ustart = ustart-addedtime;
%     uend = uend+addedtime;
%     uspec = zspec(:,ustart:uend)';
%     newspec = resample(uspec,P,Q);
    alleventspec(:,:,u) = newspec;
    alleventLFP(:,u) = newLFP;
    newspec = [1 1];
end


meaneventspec = mean(alleventspec,3);


%remove ends...
%% Figure
% figure
%         hold on
%         imagesc(t_scaled,log2(freqs),meanUP')
%         plot([0 0],log2([freqs(1) freqs(end)]),'w',...
%             [1 1],log2([freqs(1) freqs(end)]),'w')
%         set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
%         set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
%         set(gca,'XTick',[0 1])
%         axis xy
%         colorbar
%         %caxis([-1 1])
%         xlim([-0.7 1.7]);ylim(log2([freqs(1) freqs(end)]));
%         xlabel('Normalized UP time')
%         ylabel('f (Hz)')

end

