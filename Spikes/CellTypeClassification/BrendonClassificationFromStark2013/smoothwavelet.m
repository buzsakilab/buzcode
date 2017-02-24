function swave=smoothwavelet(wave,dt,period,dj,scale)
% Smoothing as in the appendix of Torrence and Webster "Inter decadal changes in the ENSO-Monsoon System" 1998
% 
% used in wavelet coherence calculations
% 
%
% Only applicable for the Morlet wavelet. 
%
% (C) Aslak Grinsted 2002-2005
%

% -------------------------------------------------------------------------
%   Copyright (C) 2002-2005, Aslak Grinsted
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.

% TODO: take mother argument
%


n=size(wave,2);

%swave=zeros(size(wave));
twave=zeros(size(wave));

% %filter in time:.... 
% for i=1:size(wave,1)
%     sc=period(i)/dt; % time/cycle / time/sample = samples/cycle
%     t=(-round(sc*3):round(sc*3))*dt;
%     f=exp(-t.^2/(2*scale(i)^2));
%     f=f/sum(f); %filter must have unit weight
%     
%     smooth=conv(wave(i,:),f); %slowest line of code in the wtcsig calculation. should be done like in wavelet with fft and ifft.
%     cutlen=(length(t)-1)*.5;
%     twave(i,:)=smooth((cutlen+1):(end-cutlen)); %remove paddings
% end
% 
%filter in time:....
% 
% qwave=twave;

%zero-pad to power of 2... Speeds up fft calcs if n is large
npad=2.^ceil(log2(n));

k = 1:fix(npad/2);
k = k.*((2.*pi)/npad);
k = [0., k, -k(fix((npad-1)/2):-1:1)];

k2=k.^2;
snorm=scale./dt;
for ii=1:size(wave,1)
    F=exp(-.5*(snorm(ii)^2)*k2); %Thanks to Bing Si for finding a bug here.
    smooth=ifft(F.*fft(wave(ii,:),npad));
    twave(ii,:)=smooth(1:n);
end

if isreal(wave)
    twave=real(twave); %-------hack-----------
end

%scale smoothing (boxcar with width of .6)

%
% TODO: optimize. Because this is done many times in the monte carlo run.
%


dj0=0.6;
dj0steps=dj0/(dj*2);
% for ii=1:size(twave,1)
%     number=0;
%     for l=1:size(twave,1);
%         if ((abs(ii-l)+.5)<=dj0steps)
%             number=number+1;
%             swave(ii,:)=swave(ii,:)+twave(l,:);
%         elseif ((abs(ii-l)+.5)<=(dj0steps+1))
%             fraction=mod(dj0steps,1);
%             number=number+fraction;
%             swave(ii,:)=swave(ii,:)+twave(l,:)*fraction;
%         end
%     end
%     swave(ii,:)=swave(ii,:)/number;
% end

kernel=[mod(dj0steps,1); ones(2 * round(dj0steps)-1,1); ...
   mod(dj0steps,1)]./(2*round(dj0steps)-1+2*mod(dj0steps,1));
swave=conv2(twave,kernel,'same'); %thanks for optimization by Uwe Graichen 

%swave=twave;