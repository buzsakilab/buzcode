function L = SpkTrainHDValuation(pf,S,x1,ang,testEp,trainingEp)

dAng= Data(Restrict(ang,testEp));
t = Range(Restrict(ang,testEp));

f = zeros(length(t),1);
r = rate(S,trainingEp);
S = Restrict(S,testEp);

xx = floor((length(x1)-1)*contrast(dAng))+1;

pf = pf(:);
try
    f = pf(xx);
catch
    keyboard
end

dt = median(diff(t));
ft = f;
ft = tsd(t,ft/r);
fs = Restrict(ft,S);
d = Data(fs);

intF = sum(f-r).*(dt/10000); %express it in sec!
d = sum(log(d(d>0)));
L = d-intF;






