function [h px py] = IcsiStat(rgx,rgy,icsi,tR,tC,tMax)


% USAGE:
%    [h px py] = IcsiStat(rgx,rgy,icsi,tR,tC,tMax)
%    
% IcsiStat computes some statistics on the common inter-spike interval
% (based on Fee et al., 1996)
% INPUTs:
%     rgx: 1st spike train
%     rgy: 2nd spike train
%     icsi: common isi comted with Icsi()
%     tR: refractory period (in seconds)
%     tC: censored time (in seconds)
%     tMax: maximal time considered for ISI (in seconds)
%     
% OUTPUT:
%     h: binary value telling of the test was successfully passed or not.
%     
% Adrien Peyrache, 2012

t = [0:0.0001:tMax];
Sx = zeros(length(t),1);
Sy = zeros(length(t),1);

isiX = diff(rgx);
isiX = isiX(isiX<tMax);
Mx = length(isiX);

isiY = diff(rgy);
isiY = isiY(isiY<tMax);
My = length(isiY);

icsi = icsi(icsi<tMax);
Mxy = length(icsi);

for ii=1:length(t)

    Fxy = sum(icsi<t(ii))/Mxy;
    Fx = sum(isiX<t(ii))/Mx;
    Fy = sum(isiY<t(ii))/My;
    if isnan(Fxy)
        Fxy=0;
    end
     if isnan(Fy)
        Fy=0;
     end
     if isnan(Fx)
        Fx=0;
    end
    Sx(ii) = sqrt(Mxy*Mx/(Mxy+Mx))*(Fxy-Fx);
    Sy(ii) = sqrt(Mxy*My/(Mxy+My))*(Fxy-Fy);    
end

f = (tR-tC)/(tMax-tC);
px = P(max(Sx(t<2)),f);
py = P(max(Sy(t<2)),f);

h = px<0.95 & py<0.95;

if 0
figure(1),clf
plot(Sx)
hold on
plot(Sy,'r')
end

function P=P(t,f)
P = 0.5*(1+erf(t/(sqrt(2*f*(1-f)))))-0.5*exp(-2*t.^2).*(1-erf(2*(1-2*f)*t/sqrt(2*f*(1-f))));
