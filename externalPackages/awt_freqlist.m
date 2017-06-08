function [wt,freqlist,psi_array] = awt_freqlist(x,Fs,freqlist,type,xi)
%   awt_freqlist  analytical wavelet transform, where one can specify the list of desired frequencies 
%   
%   [wt,freqlist,psi_array] = awt_freqlist(x,Fs,freqlist,type,xi)
%
%   Inputs:
%       x           the signal to be analyzed
%       Fs          the sampling frequency
%       freqlist    list of frequencies at which to compute wt (optional)
%                   (or set to [] for automatic definition)
%       type        type of wavelet to use (Gabor, Lusin or Sombrero)
%       xi          the number of oscillations parameter
%   Outputs:
%       wt: time-frequency image
%       freqlist: useful in case of automatic definition of frequencies
%       psi_array : array of analysis functions (complex values)
%
%  Maureen Clerc, Christian Benar, october 2007
%  modified from awt from wavelab toolbox

% History of changes
% 1/11/2007: (cgb) psi_array: output in complex 
% 3/06/2008: (cgb) init of psi_array to size of wt

n = length(x);
sigma2 = 1;
x = x(:);
omega = [(0:n/2) (-ceil(n/2)+1:-1)].*Fs/n; % CGB added ceil: TO BE CHECKED
omega = omega(:);

fftx = fft(x);

% to compute min and max center frequency values:
%tolerance = 1.5; % such that an acceptable "zero" is exp(-tolerance^2/2)
tolerance = 0.5; % cgb
%tolerance = 1; % cgb

if nargin<2 
    Fs = 1;
end
if nargin<4
    type = 'Gabor';
end
if nargin<5
    xi = 5; % only useful for Gabor
end

% compute min and max valid center frequencies, according to wavelet type
switch(type)
    case 'Gabor'
        mincenterfreq = 2*tolerance*sqrt(sigma2)*Fs*xi./n;
        maxcenterfreq = Fs*xi/(xi+tolerance/sqrt(sigma2));
    case 'Lusin1'
        % to be computed
    case 'Sombrero'
        [mincenterfreq,maxcenterfreq] = FindSombreroBounds(tolerance,sigma2,n,Fs);
end
if  (nargin<3 || isempty(freqlist) )
    nvoice = 12;
    freqlist= 2.^(log2(mincenterfreq):1/nvoice:log2(maxcenterfreq));
end

switch(type)
    case 'Gabor'
        s_array = xi./freqlist;
        minscale = xi./maxcenterfreq;
        maxscale = xi./mincenterfreq;
    case 'Lusin1'
        % to be computed
    case 'Sombrero'
        s_array = sqrt(2/sigma2)./freqlist;
        minscale = sqrt(2/sigma2)/maxcenterfreq;
        maxscale = sqrt(2/sigma2)/mincenterfreq;
end
nscale = length(freqlist);
wt = zeros(n,nscale);
scaleindices=find(s_array(:)'>=minscale & s_array(:)'<=maxscale);
%psi_array=zeros(n,length(scaleindices));
psi_array=zeros(n,nscale);
for kscale=scaleindices
    s=s_array(kscale);
    switch(type)
        case 'Gabor'
            freq =  (s .* omega  - xi);
            Psi = realpow(4.*pi.*sigma2,1/4)*sqrt(s) *exp(-sigma2/2*freq.*freq);
            
        %    Psi = Psi .* (omega>0);
        case 'Lusin1'
            % to be completed
        case 'Sombrero'
            freq = s.*omega;
            Psi = -realpow(pi*sigma2^5,1/4)*sqrt(8/3)*freq.^2.*sqrt(s).*exp(-sigma2/2*freq.*freq);
            Psi = Psi .* (omega>0);
    end;

    wt(1:n,kscale) = ifft(fftx.*Psi);
    psi_array(:,kscale)=ifft(Psi);
end


function [mincenterfreq,maxcenterfreq] = FindSombreroBounds(tol,sigma2,n,Fs)
s = 0.02; % an intermediate scale according to which the min and max scale will be computed
t = (1:n)/Fs./s;
psi = 2*realpow(pi,-1/4)/sqrt(3).*(t.^2-1).*exp(-t.^2./(2*sigma2));
psi = psi./max(psi);
% in the decreasing part of psi, find first sub-threshold
i = find(abs(psi)<exp(-tol^2/2) & psi-circshift(psi,[0 1])<0,1,'first');
smax = s*n/i; % max scale at which sombrero remains within time window
mincenterfreq = sqrt(2/sigma2)/smax;
omega = s*(0:n/2).*Fs/n;
hatpsi = realpow(sigma2*pi,1/4)*sqrt(8/3)*sigma2*omega.^2.*exp(-omega.^2.*sigma2./2);
hatpsi = hatpsi./max(hatpsi);
i = find(hatpsi>exp(-tol^2/2) & hatpsi-circshift(hatpsi,[0 1])<0,1,'last');
smin = i/(n/2)*s; %l'echelle min a laquelle le sombrero ne sort pas de la fenetre frequentielle
maxcenterfreq = sqrt(2/sigma2)/smin;
