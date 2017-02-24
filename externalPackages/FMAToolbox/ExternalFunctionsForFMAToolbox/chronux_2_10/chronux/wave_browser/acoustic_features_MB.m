function [Feat,S,t,f]=acoustic_features_MB(data,movingwin,params)

%   Usage:
%   [Feat,t]=acoustic_features_MB(data,movingwin,params);
%   Input: 
%   Note units have to be consistent. Thus, if movingwin is in seconds, Fs
%   has to be in Hz. see chronux.m for more information.
%         data        Time series -- required
%         movingwin         (in the form [window winstep] i.e length of moving
%                                                   window and step size)
%                                                   Note that units here have
%                                                   to be consistent with
%                                                   units of Fs - required
%         params: structure with fields tapers, pad, Fs, fpass
%         - optional
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
%
%  	        pad		    (padding factor for the FFT) - optional. Defaults to 0.  
%  			      	 e.g. For N = 500, if PAD = 0, we pad the FFT 
%  			      	 to 512 points; if PAD = 2, we pad the FFT
%  			      	 to 2048 points, etc.
%             Fs   (sampling frequency) - optional. Default 1.
%             fpass    (frequency band to be used in the calculation in the form
%                                     [fmin fmax])- optional. 
%                                     Default all frequencies between 0 and Fs/2
% 
%   Output:
%         Featt   Features: 3-dim time series, <S>, <log(S)> and <f S>/<S>.
%         first two averages computed over fpass.
%         t       (times)

params1=params;
fpass=params.fpass;fpass1=fpass;
fpass1(1)=0;
params1.fpass=fpass1;
[S,t,f]=mtspecgramc(diff(data),movingwin,params1);
Feat=zeros(length(t),3);
pass=floor(fpass/params.Fs*length(f))+1;
Feat(:,1)=mean(S(:,pass(1):pass(2)),2);
Feat(:,2)=mean(log(S(:,pass(1):pass(2))),2);
f=f(:)';
freq=repmat(f,length(t),1);
% Feat(:,3)=mean(freq.*S,2)./mean(S,2);
Feat(:,3)=max(S,[],2)./median(S,2);
