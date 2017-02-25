function [S,t,f,Serr]=mtspecgramc(data,movingwin,params)
% Multi-taper time-frequency spectrum - continuous process
%
% Usage:
% [S,t,f,Serr]=mtspecgramc(data,movingwin,params)
% Input: 
% Note units have to be consistent. Thus, if movingwin is in seconds, Fs
% has to be in Hz. see chronux.m for more information.
%       data        (in form samples x channels/trials) -- required
%       movingwin         (in the form [window winstep] i.e length of moving
%                                                 window and step size)
%                                                 Note that units here have
%                                                 to be consistent with
%                                                 units of Fs - required
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%       - optional
%           tapers : precalculated tapers from dpss or in the one of the following
%                    forms: 
%                   (1) A numeric vector [TW K] where TW is the
%                       time-bandwidth product and K is the number of
%                       tapers to be used (less than or equal to
%                       2TW-1). 
%                   (2) A numeric vector [W T p] where W is the
%                       bandwidth, T is the duration of the data and p 
%                       is an integer such that 2TW-p tapers are used. In
%                       this form there is no default i.e. to specify
%                       the bandwidth, you have to specify T and p as
%                       well. Note that the units of W and T have to be
%                       consistent: if W is in Hz, T must be in seconds
%                       and vice versa. Note that these units must also
%                       be consistent with the units of params.Fs: W can
%                       be in Hz if and only if params.Fs is in Hz.
%                       The default is to use form 1 with TW=3 and K=5
%                   Note that T has to be equal to movingwin(1).
%
%	        pad		    (padding factor for the FFT) - optional. Defaults to 0.  
%			      	 e.g. For N = 500, if PAD = 0, we pad the FFT 
%			      	 to 512 points; if PAD = 2, we pad the FFT
%			      	 to 2048 points, etc.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%           trialave (average over trials/channels when 1, don't average when 0) - optional. Default 0
% Output:
%       S       (spectrum in form time x frequency x channels/trials if trialave=0; in the form time x frequency if trialave=1)
%       t       (times)
%       f       (frequencies)
%       Serr    (error bars) only for err(1)>=1

if nargin < 2; error('Need data and window parameters'); end;
if nargin < 3; params=[]; end;

if length(params.tapers)==3 & movingwin(1)~=params.tapers(2);
    error('Duration of data in params.tapers is inconsistent with movingwin(1), modify params.tapers(2) to proceed')
end

[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if nargout > 3 && err(1)==0; 
%   Cannot compute error bars with err(1)=0. change params and run again.
    error('When Serr is desired, err(1) has to be non-zero.');
end;

N=size(data,1);
Nwin=round(Fs*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*Fs); % number of samples to step through
nfft=2^(nextpow2(Nwin)+pad);

[f,findx]=getfgrid(Fs,nfft,fpass);
tapers=dpsschk(tapers,Nwin,Fs); % check tapers
[NC C]=size(data); % size of data
[NK K]=size(tapers); % size of tapers

tapers=tapers(:,:,ones(1,C)); % add channel indices to tapers

if NK~=Nwin ; error('length of tapers is incompatible with length of data'); end;

winstart=1:Nstep:N-Nwin+1;
nw=length(winstart); 

S = zeros(nw,length(f),C);

for n=1:nw;
    J=fft((permute(data(winstart(n):winstart(n)+Nwin-1,:,ones(1,K)),[1 3 2]) .* tapers)  ,nfft)/Fs;   % fft of projected data
    J=J(findx,:,:);
    S(n,:,:)=squeeze(mean(conj(J).*J,2));
end

if nargout==4;Serr=squeeze(Serr);end;
winmid=winstart+round(Nwin/2);
t=winmid/Fs;
