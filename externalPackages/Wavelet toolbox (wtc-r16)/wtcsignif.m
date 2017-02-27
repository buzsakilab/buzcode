function wtcsig=wtcsignif(mccount,ar1,dt,n,pad,dj,s0,j1,mother,cutoff)
% Wavelet Coherence Significance Calculation (Monte Carlo)
%
% wtcsig=wtcsignif(mccount,ar1,dt,n,pad,dj,s0,j1,mother,cutoff)
%
% mccount: number of time series generations in the monte carlo run (the greater the better)
% ar1: a vector containing the 2 AR(1) coefficients. example: ar1=[.7 .6].
% dt,pad,dj,s0,j1,mother: see wavelet help... 
% n: length of each generated timeseries. (obsolete) 
%
% cutoff: (obsolete)
%
% RETURNED
% wtcsig: the 95% significance level as a function of scale... (scale,sig95level)
%
% (C) Aslak Grinsted 2002-2005
%

% -------------------------------------------------------------------------
%   Copyright (C) 2002-2005, Aslak Grinsted
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.
persistent mypath
if isempty(mypath)
    mypath=strrep(which('wtcsignif'),'wtcsignif.m','');
end




% TODO: also make XWT Signif calculation since it is cheap. (?)

%we don't need to do the monte carlo if we have a cached
%siglevel for ar1s that are almost the same. (see fig4 in Grinsted et al., 2004)
aa=round(atanh(ar1(:)')*4); %this function increases the sensitivity near 1 & -1
aa=abs(aa)+.5*(aa<0); %only positive numbers are allowed in the checkvalues (because of log)

% do a check that it is not the same as last time... (for optimization purposes)
checkvalues=[aa dj s0/dt j1 double(mother)]; %n & pad are not important.
%also the resolution is not important.

checkhash=['' mod(sum(log(checkvalues+1)*127),25)+'a' mod(sum(log(checkvalues+1)*54321),25)+'a'];

cachefilename=[mypath 'wtcsignif-cached-' checkhash '.bnm'];

%the hash is used to distinguish cache files.
try
    [lastmccount,lastcheckvalues,lastwtcsig]=loadbnm(cachefilename);
    if (lastmccount>=mccount)&(isequal(single(checkvalues),lastcheckvalues)) %single is important because bnm is single precision.
        wtcsig=lastwtcsig;
        return
    end
catch
end

%choose a n so that largest scale have atleast some part outside the coi
ms=s0*(2^(j1*dj))/dt; %maxscale in units of samples
n=ceil(ms*6);

warned=0;
%precalculate stuff that's constant outside the loop
%d1=ar1noise(n,1,ar1(1),1);
d1=rednoise(n,ar1(1),1);
[W1,period,scale,coi] = wavelet(d1,dt,pad,dj,s0,j1,mother);
outsidecoi=zeros(size(W1));
for s=1:length(scale)
    outsidecoi(s,:)=(period(s)<=coi);
end
sinv=1./(scale');
sinv=sinv(:,ones(1,size(W1,2)));

if mccount<1
    wtcsig=scale';
    wtcsig(:,2)=.71; %pretty good 
    return
end

sig95=zeros(size(scale));

maxscale=1;
for s=1:length(scale)
    if any(outsidecoi(s,:)>0)
        maxscale=s;
    else
        sig95(s)=NaN;
        if ~warned
            warning('Long wavelengths completely influenced by COI. (suggestion: set n higher, or j1 lower)'); warned=1;
        end
    end
end

%PAR1=1./ar1spectrum(ar1(1),period');
%PAR1=PAR1(:,ones(1,size(W1,2)));
%PAR2=1./ar1spectrum(ar1(2),period');
%PAR2=PAR2(:,ones(1,size(W1,2)));

nbins=1000;
wlc=zeros(length(scale),nbins);

wbh = waitbar(0,['Running Monte Carlo (significance)... (H=' checkhash ')'],'Name','Monte Carlo (WTC)');
for ii=1:mccount
    waitbar(ii/mccount,wbh);
%     d1=ar1noise(n,1,ar1(1),1);    
%     d2=ar1noise(n,1,ar1(2),1);    
    d1=rednoise(n,ar1(1),1);    
    d2=rednoise(n,ar1(2),1);    
    [W1,period,scale,coi] = wavelet(d1,dt,pad,dj,s0,j1,mother);
    [W2,period,scale,coi] = wavelet(d2,dt,pad,dj,s0,j1,mother);
%    W1=W1.*PAR1; %whiten spectra
%    W2=W2.*PAR2;
    sWxy=smoothwavelet(sinv.*(W1.*conj(W2)),dt,period,dj,scale);
    Rsq=abs(sWxy).^2./(smoothwavelet(sinv.*(abs(W1).^2),dt,period,dj,scale).*smoothwavelet(sinv.*(abs(W2).^2),dt,period,dj,scale));
    
    for s=1:maxscale
        cd=Rsq(s,find(outsidecoi(s,:)));
        cd=max(min(cd,1),0);
        cd=floor(cd*(nbins-1))+1;
        for jj=1:length(cd)
            wlc(s,cd(jj))=wlc(s,cd(jj))+1;
        end
    end
end
close(wbh);



for s=1:maxscale
    rsqy=((1:nbins)-.5)/nbins;
    ptile=wlc(s,:);
    idx=find(ptile~=0);
    ptile=ptile(idx);
    rsqy=rsqy(idx);
    ptile=cumsum(ptile);
    ptile=(ptile-.5)/ptile(end);
    sig95(s)=interp1(ptile,rsqy,.95);
end
wtcsig=[scale' sig95'];

if any(isnan(sig95))&(~warned)
    warning(sprintf('Sig95 calculation failed. (Some NaNs)'))
else
    savebnm(cachefilename,mccount,checkvalues,wtcsig); %save to a cache....
end


