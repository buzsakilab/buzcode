function [S,t,f,R,Serr]=mtspecgrampt(data,movingwin,params,fscorr)
% Multi-taper time-frequency spectrum - point process times
%
% Usage:
%
% [S,t,f,R,Serr]=mtspecgrampt(data,movingwin,params,fscorr)
% Input: 
%       data        (structure array of spike times with dimension channels/trials; 
%                   also accepts 1d array of spike times) -- required
%       movingwin         (in the form [window,winstep] i.e length of moving
%                                                 window and step size.
%                                                 
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%       - optional
%           tapers : precalculated tapers from dpss or in the one of the following
%                    forms: 
%                    (1) A numeric vector [TW K] where TW is the
%                        time-bandwidth product and K is the number of
%                        tapers to be used (less than or equal to
%                        2TW-1). 
%                    (2) A numeric vector [W T p] where W is the
%                        bandwidth, T is the duration of the data and p 
%                        is an integer such that 2TW-p tapers are used. In
%                        this form there is no default i.e. to specify
%                        the bandwidth, you have to specify T and p as
%                        well. Note that the units of W and T have to be
%                        consistent: if W is in Hz, T must be in seconds
%                        and vice versa. Note that these units must also
%                        be consistent with the units of params.Fs: W can
%                        be in Hz if and only if params.Fs is in Hz.
%                        The default is to use form 1 with TW=3 and K=5
%                     Note that T has to be equal to movingwin(1).
%
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%           trialave (average over trials/channels when 1, don't average when 0) - optional. Default 0
%       fscorr   (finite size corrections, 0 (don't use finite size corrections) or 
%                1 (use finite size corrections) - optional
%                (available only for spikes). Defaults 0.
%
% Output:
%       S       (spectrogram with dimensions time x frequency x channels/trials if trialave=0; 
%               dimensions time x frequency if trialave=1)
%       t       (times)
%       f       (frequencies)
%
%       Serr    (error bars) - only if err(1)>=1
%
% This is an optimized version of Chronux 2.11, which pre-calculates tapers and their fft
% before looping through the spectrogram windows. Optimization by Ralitsa Todorova 2016.

if nargin < 2; error('Need data and window parameters'); end;
if nargin < 3; params=[]; end;

[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if length(params.tapers)==3 & movingwin(1)~=params.tapers(2);
    error('Duration of data in params.tapers is inconsistent with movingwin(1), modify params.tapers(2) to proceed')
end

data=change_row_to_column(data);
if isstruct(data); Ch=length(data); end;
if nargin < 4 || isempty(fscorr); fscorr=0; end;
if nargout > 4 && err(1)==0; error('Cannot compute errors with err(1)=0'); end;

[mintime,maxtime]=minmaxsptimes(data);
tn=(mintime+movingwin(1)/2:movingwin(2):maxtime-movingwin(1)/2);
Nwin=round(Fs*movingwin(1)); % number of samples in window
nfft=max(2^(nextpow2(Nwin)+pad),Nwin);
f=getfgrid(Fs,nfft,fpass); Nf=length(f);
params.tapers=dpsschk(tapers,Nwin,Fs); % check tapers
nw=length(tn);

if trialave
    S = zeros(nw,Nf);
    R = zeros(nw,1);
    if nargout==4; Serr=zeros(2,nw,Nf); end;
else
    S = zeros(nw,Nf,Ch);
    R = zeros(nw,Ch);
    if nargout==4; Serr=zeros(2,nw,Nf,Ch); end;
end

tapers=dpsschk(tapers,Nwin,Fs); % check tapers
[f,findx]=getfgrid(Fs,nfft,fpass); % get frequency grid for evaluation
K=size(tapers,2); % number of tapers
nfreq=length(f); % number of frequencies
H=fft(tapers,nfft,1);  % fft of tapers
H=H(findx,:); % restrict fft of tapers to required frequencies
w=2*pi*f; % angular frequencies at which ft is to be evaluated

for n=1:nw;
    t=linspace(tn(n)-movingwin(1)/2,tn(n)+movingwin(1)/2,Nwin);
    datawin=extractdatapt(data,[t(1) t(end)]);
    datawin=change_row_to_column(datawin);
    if isstruct(datawin); C=length(datawin); else C=1; end% number of channels
    Nsp=zeros(1,C); Msp=zeros(1,C);
    for ch=1:C;
        if isstruct(datawin);
            fnames=fieldnames(datawin);
            eval(['dtmp=datawin(ch).' fnames{1} ';'])
            indx=find(dtmp>=min(t)&dtmp<=max(t));
            if ~isempty(indx); dtmp=dtmp(indx);
            end;
        else
            dtmp=datawin;
            indx=find(dtmp>=min(t)&dtmp<=max(t));
            if ~isempty(indx); dtmp=dtmp(indx);
            end;
        end;
    end
    Nsp(ch)=length(dtmp);
    Msp(ch)=Nsp(ch)/length(t);
    if Msp(ch)~=0;
        data_proj=interp1(t',tapers,dtmp);
        exponential=exp(-i*w'*(dtmp-t(1))');
        J(:,:,ch)=exponential*data_proj-H*Msp(ch);
    else
        J(1:nfreq,1:K,ch)=0;
    end;
    s=squeeze(mean(conj(J).*J,2));
    if trialave; s=squeeze(mean(s,2));Msp=mean(Msp);end;
    r=Msp*Fs;
    if nargout==5;
        if fscorr==1;
            serr=specerr(S,J,err,trialave,Nsp);
        else
            serr=specerr(S,J,err,trialave);
        end
        Serr(1,n,:,:)=squeeze(serr(1,:,:));
        Serr(2,n,:,:)=squeeze(serr(2,:,:));
    end;
    S(n,:,:)=s;
    R(n,:)=r;
end;
t=tn;
S=squeeze(S); R=squeeze(R); if nargout==5; Serr=squeeze(Serr);end


