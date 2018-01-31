
function [wave,period,scale,coi] = ...
    	wavelet_mod(Y,dt,pad,frecs,param)

% Calcula la transformada wavelet MORLET
% entrada
% Y= senal de entrada nptosx1
% dt = 1/Fs
% pad = 1 si queremos padding con 0s hasta que la longitud de la señal sea potencia de 2
% frecs = vector con las frecuencias para las que queremos calcular wavlet
% param =parametro morlet (defecto =6)

n1 = length(Y);
mother='MORLET';
if (nargin<5), param = 6.; end


%....construct time series to analyze, pad if necessary
x(1:n1) = Y - mean(Y);
if (pad == 1)
	base2 = fix(log(n1)/log(2) + 0.4999);   % power of 2 nearest to N
	x = [x,zeros(1,2^(base2+1)-n1)];
end
n = length(x);

%....construct wavenumber array used in transform [Eqn(5)]
k = 1:fix(n/2);
k = k.*((2.*pi)/(n*dt));
k = [0., k, -k(fix((n-1)/2):-1:1)];

%....compute FFT of the (padded) time series
f = fft(x);    % [Eqn(3)]

%....construct SCALE array & empty PERIOD & WAVE arrays
fourier_factor = (4*pi)/(param + sqrt(2 + param^2)); % Scale-->Fourier [Sec.3h]
scale=1./(fourier_factor*frecs);
J1=length(scale)-1;

wave = zeros(J1+1,n);  % define the wavelet array
% loop through all scales and compute transform
for a1 = 1:J1+1
	[daughter,fourier_factor,coi,~]=wave_bases(mother,k,scale(a1),param);	
	wave(a1,:) = ifft(f.*daughter);  % wavelet transform[Eqn(4)]
end

period = fourier_factor*scale;
coi = coi*dt*[1E-5,1:((n1+1)/2-1),fliplr((1:(n1/2-1))),1E-5];  % COI [Sec.3g]
wave = wave(:,1:n1);  % get rid of padding before returning

return
